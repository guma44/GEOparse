"""
Classes that represent different GEO entities
"""

import os
import re
import abc
import gzip
import glob
import json
import time
import platform
import subprocess
import numpy as np
from pandas import DataFrame, concat

try:
    from urllib.error import HTTPError
except ImportError:
    from urllib2 import HTTPError
from six import iteritems, itervalues

from . import utils
from .logger import logger


class DataIncompatibilityException(Exception):
    pass


class NoMetadataException(Exception):
    pass


class NoSRARelationException(Exception):
    pass


class NoSRAToolkitException(Exception):
    pass


class BaseGEO(object):
    __metaclass__ = abc.ABCMeta

    geotype = None

    def __init__(self, name, metadata):
        """Initialize base GEO object.

        Args:
            name (:obj:`str`): Name of the object.
            metadata (:obj:`dict`): Metadata information.

        Raises:
            TypeError: Metadata should be a dict.
        """

        if not isinstance(metadata, dict):
            raise TypeError("Metadata should be a dictionary not a %s" % str(
                type(metadata)))

        self.name = name
        self.metadata = metadata
        self.relations = {}
        if 'relation' in self.metadata:
            for relation in self.metadata['relation']:
                tmp = re.split(r':\s+', relation)
                relname = tmp[0]
                relval = tmp[1]

                if relname in self.relations:
                    self.relations[relname].append(relval)
                else:
                    self.relations[relname] = [relval]

    def get_metadata_attribute(self, metaname):
        """Get the metadata attribute by the name.

        Args:
            metaname (:obj:`str`): Name of the attribute

        Returns:
            :obj:`list` or :obj:`str`: Value(s) of the requested metadata
                attribute

        Raises:
            NoMetadataException: Attribute error
            TypeError: Metadata should be a list
        """
        metadata_value = self.metadata.get(metaname, None)
        if metadata_value is None:
            raise NoMetadataException(
                "No metadata attribute named %s" % metaname)
        if not isinstance(metadata_value, list):
            raise TypeError("Metadata is not a list and it should be.")

        if len(metadata_value) > 1:
            return metadata_value
        else:
            return metadata_value[0]

    def get_accession(self):
        """Return accession ID of the sample.

        Returns:
            :obj:`str`: GEO accession ID
        """
        return self.get_metadata_attribute("geo_accession")

    def get_type(self):
        """Get the type of the GEO.

        Returns:
            :obj:`str`: Type attribute of the GEO
        """
        try:
            return self.get_metadata_attribute("type")
        except NoMetadataException:
            return None

    def _get_metadata_as_string(self):
        """Get the metadata as SOFT formatted string."""
        metalist = []
        for metaname, meta in iteritems(self.metadata):
            message = "Single value in metadata dictionary should be a list!"
            assert isinstance(meta, list), message
            for data in meta:
                if data:
                    metalist.append("!%s_%s = %s" % (self.geotype.capitalize(),
                                                     metaname, data))
        return "\n".join(metalist)

    def show_metadata(self):
        """Print metadata in SOFT format."""
        print(self._get_metadata_as_string())

    def to_soft(self, path_or_handle, as_gzip=False):
        """Save the object in a SOFT format.

        Args:
            path_or_handle (:obj:`str` or :obj:`file`): Path or handle to
                output file
            as_gzip (:obj:`bool`): Save as gzip
        """
        if isinstance(path_or_handle, str):
            if as_gzip:
                with gzip.open(path_or_handle, 'wt') as outfile:
                    outfile.write(self._get_object_as_soft())
            else:
                with open(path_or_handle, 'w') as outfile:
                    outfile.write(self._get_object_as_soft())
        else:
            path_or_handle.write(self._get_object_as_soft())

    @abc.abstractmethod
    def _get_object_as_soft(self):
        """Get the object as SOFT formatted string."""
        raise NotImplementedError("Method not implemented")

    def __str__(self):
        return str("<%s: %s>" % (self.geotype, self.name))

    def __repr__(self):
        return str("<%s: %s>" % (self.geotype, self.name))


class SimpleGEO(BaseGEO):
    __metaclass__ = abc.ABCMeta

    def __init__(self, name, metadata, table, columns):
        """Initialize simple GEO object.

        Args:
            name (:obj:`str`): Name of the object
            metadata (:obj:`dict`): Metadata information
            table (:obj:`pandas.DataFrame`): Table with the data from SOFT file
            columns (:obj:`pandas.DataFrame`): Description of the columns,
                number of columns, order and names represented as index
                in this DataFrame has to be the same as table.columns.

        Raises:
            :obj:`ValueError`: Table should be a DataFrame
            :obj:`ValueError`: Columns' description should be a DataFrame
            :obj:`DataIncompatibilityException`: Columns are wrong
            :obj:`ValueError`: Description has to be present in columns
        """
        if not isinstance(table, DataFrame):
            raise ValueError(("Table data should be an instance of "
                              "pandas.DataFrame not %s") % str(type(table)))
        if not isinstance(columns, DataFrame):
            raise ValueError(("Columns description should be an instance of "
                              "pandas.DataFrame not %s") % str(type(columns)))

        BaseGEO.__init__(self, name=name, metadata=metadata)

        self.table = table
        self.columns = columns
        columns_are_correct = False
        if self.columns.index.tolist() != self.table.columns.tolist():
            if not self.columns.index.is_unique:
                # try to correct duplicate index in the same way pandas is
                # doing the columns
                logger.warning(
                    "Detected duplicated columns in %s %s. Correcting.\n" % (
                        self.geotype,
                        self.name))
                indices = {}
                new_index = []
                for idx in self.columns.index:
                    if idx not in indices:
                        indices[idx] = 0
                        new_index.append(idx)
                    else:
                        indices[idx] += 1
                        new_index.append("%s.%i" % (idx, indices[idx]))
                self.columns.index = new_index
                if self.columns.index.tolist() == self.table.columns.tolist():
                    columns_are_correct = True
            if not columns_are_correct:
                # if the columns are still not correct check the order.
                if sorted(self.columns.index.tolist()) == sorted(
                    self.table.columns.tolist()):
                    logger.warning(
                        "Data columns in %s %s are not in order. Reordering.\n" % (
                            self.geotype,
                            self.name))
                    self.columns = self.columns.ix[self.table.columns]
                    if self.columns.index.tolist() == self.table.columns.tolist():
                        columns_are_correct = True
        else:
            columns_are_correct = True

        if not columns_are_correct:
            rows_in_columns = ", ".join(self.columns.index.tolist())
            columns_in_table = ", ".join(self.table.columns.tolist())
            raise DataIncompatibilityException(
                "\nData columns do not match columns description index in %s\n" % (
                    self.name) +
                "Columns in table are: %s\n" % columns_in_table +
                "Index in columns are: %s\n" % rows_in_columns)
        if self.columns.columns[0] != 'description':
            raise ValueError(("Columns table must contain a column named"
                              "'description'. Here columns are: %s") % (
                                 ", ".join(map(str, self.columns.columns))))

    def head(self):
        """Print short description of the object."""
        summary = list()
        summary.append("%s %s" % (self.geotype, self.name) + "\n")
        summary.append(" - Metadata:" + "\n")
        summary.append(
            "\n".join(self._get_metadata_as_string().split("\n")[:5]) + "\n")
        summary.append("\n")
        summary.append(" - Columns:" + "\n")
        summary.append(self.columns.to_string() + "\n")
        summary.append("\n")
        summary.append(" - Table:" + "\n")
        summary.append(
            "\t".join(["Index"] + self.table.columns.tolist()) + "\n")
        summary.append(self.table.head().to_string(header=None) + "\n")
        summary.append(" " * 40 + "..." + " " * 40 + "\n")
        summary.append(" " * 40 + "..." + " " * 40 + "\n")
        summary.append(" " * 40 + "..." + " " * 40 + "\n")
        summary.append(self.table.tail().to_string(header=None) + "\n")
        return "\n".join(summary)

    def show_columns(self):
        """Print columns in SOFT format."""
        print(self.columns)

    def show_table(self, number_of_lines=5):
        """Show few lines of the table the table as pandas.DataFrame.

        Args:
            number_of_lines (:obj:`int`): Number of lines to show. Defaults to 5.
        """
        print(self.table.head(number_of_lines))

    def _get_object_as_soft(self):
        """Get the object as SOFT formated string."""
        soft = ["^%s = %s" % (self.geotype, self.name),
                self._get_metadata_as_string(),
                self._get_columns_as_string(),
                self._get_table_as_string()]
        return "\n".join(soft)

    def _get_table_as_string(self):
        """Get table as SOFT formated string."""
        tablelist = []
        tablelist.append("!%s_table_begin" % self.geotype.lower())
        tablelist.append("\t".join(self.table.columns))
        for idx, row in self.table.iterrows():
            tablelist.append("\t".join(map(str, row)))
        tablelist.append("!%s_table_end" % self.geotype.lower())
        return "\n".join(tablelist)

    def _get_columns_as_string(self):
        """Returns columns  as SOFT formated string."""
        columnslist = []
        for rowidx, row in self.columns.iterrows():
            columnslist.append("#%s = %s" % (rowidx, row.description))
        return "\n".join(columnslist)


class GSM(SimpleGEO):
    """Class that represents sample from GEO database."""

    geotype = 'SAMPLE'

    def annotate(self, gpl, annotation_column, gpl_on="ID", gsm_on="ID_REF",
                 in_place=False):
        """Annotate GSM with provided GPL

        Args:
            gpl (:obj:`pandas.DataFrame`): A Platform or DataFrame to annotate with
            annotation_column (str`): Column in a table for annotation
            gpl_on (:obj:`str`): Use this column in GSM to merge. Defaults to "ID".
            gsm_on (:obj:`str`): Use this column in GPL to merge.
                Defaults to "ID_REF".
            in_place (:obj:`bool`): Substitute table in GSM by new annotated
                table. Defaults to False.

        Returns:
            :obj:`pandas.DataFrame` or :obj:`None`: Annotated table or None


        Raises:
            TypeError: GPL should be GPL or pandas.DataFrame
        """
        if isinstance(gpl, GPL):
            annotation_table = gpl.table
        elif isinstance(gpl, DataFrame):
            annotation_table = gpl
        else:
            raise TypeError("gpl should be a GPL object or a pandas.DataFrame")

        # annotate by merging
        annotated = self.table.merge(
            annotation_table[[gpl_on, annotation_column]], left_on=gsm_on,
            right_on=gpl_on)
        del annotated[gpl_on]
        if in_place:
            self.table = annotated
            return None
        else:
            return annotated

    def annotate_and_average(self, gpl, expression_column, group_by_column,
                             rename=True,
                             force=False, merge_on_column=None, gsm_on=None,
                             gpl_on=None):
        """Annotate GSM table with provided GPL.

        Args:
            gpl (:obj:`GEOTypes.GPL`): Platform for annotations
            expression_column (:obj:`str`): Column name which "expressions"
                are represented
            group_by_column (:obj:`str`): The data will be grouped and averaged
                over this column and only this column will be kept
            rename (:obj:`bool`): Rename output column to the
                self.name. Defaults to True.
            force (:obj:`bool`): If the name of the GPL does not match the platform
                name in GSM proceed anyway. Defaults to False.
            merge_on_column (:obj:`str`): Column to merge the data
                on. Defaults to None.
            gsm_on (:obj:`str`): In the case columns to merge are different in GSM
                and GPL use this column in GSM. Defaults to None.
            gpl_on (:obj:`str`): In the case columns to merge are different in GSM
                and GPL use this column in GPL. Defaults to None.

        Returns:
            :obj:`pandas.DataFrame`: Annotated data
        """
        if gpl.name != self.metadata['platform_id'][0] and not force:
            raise KeyError("Platforms from GSM (%s) and from GPL (%s)" % (
                gpl.name, self.metadata['platform_id']) +
                           " are incompatible. Use force=True to use this GPL.")
        if merge_on_column is None and gpl_on is None and gsm_on is None:
            raise Exception("You have to provide one of the two: "
                            "merge_on_column or gpl_on and gsm_on parameters")
        if merge_on_column:
            logger.info("merge_on_column is not None. Using this option.")
            tmp_data = self.table.merge(gpl.table, on=merge_on_column,
                                        how='outer')
            tmp_data = tmp_data.groupby(group_by_column).mean()[
                [expression_column]]
        else:
            if gpl_on is None or gsm_on is None:
                raise Exception("Please provide both gpl_on and gsm_on or "
                                "provide merge_on_column only")
            tmp_data = self.table.merge(gpl.table, left_on=gsm_on,
                                        right_on=gpl_on, how='outer')
            tmp_data = tmp_data.groupby(group_by_column).mean()[
                [expression_column]]
        if rename:
            tmp_data.columns = [self.name]
        return tmp_data

    def download_supplementary_files(self, directory="./", download_sra=True,
                                     email=None, sra_kwargs=None):
        """Download all supplementary data available for the sample.

        Args:
            directory (:obj:`str`): Directory to download the data (in this directory
                function will create new directory with the files).
                Defaults to "./".
            download_sra (:obj:`bool`): Indicates whether to download SRA raw
                data too. Defaults to True.
            email (:obj:`str`): E-mail that will be provided to the Entrez.
                It is mandatory if download_sra=True. Defaults to None.
            sra_kwargs (:obj:`dict`, optional): Kwargs passed to the
                download_SRA method. Defaults to None.

        Returns:
            :obj:`dict`: A key-value pair of name and paths downloaded
        """
        directory_path = os.path.abspath(
            os.path.join(directory, "%s_%s_%s" % (
                'Supp',
                self.get_accession(),
                # the directory name cannot contain many of the signs
                re.sub(r'[\s\*\?\(\),\.;]', '_', self.metadata['title'][0]))))

        utils.mkdir_p(os.path.abspath(directory_path))
        downloaded_paths = dict()
        if sra_kwargs is None:
            sra_kwargs = {}
        for metakey, metavalue in iteritems(self.metadata):
            if 'supplementary_file' in metakey:
                assert len(metavalue) == 1 and metavalue != ''
                # logger.info("Downloading %s\n" % metavalue)
                if 'sra' in metavalue[0] and download_sra:
                    downloaded_files = self.download_SRA(email,
                                                         directory=directory,
                                                         **sra_kwargs)
                    downloaded_paths[metavalue[0]] = downloaded_files
                else:
                    download_path = os.path.abspath(os.path.join(
                        directory,
                        os.path.join(directory_path,
                                     metavalue[0].split("/")[-1])))
                    utils.download_from_url(metavalue[0], download_path)
                    downloaded_paths[metavalue[0]] = download_path
        return downloaded_paths

    def download_SRA(self, email, directory='./', **kwargs):
        """Download RAW data as SRA file.

        The files will be downloaded to the sample directory created ad hoc
        or the directory specified by the parameter. The sample has to come
        from sequencing eg. mRNA-seq, CLIP etc.

        An important parameter is a download_type. By default an SRA
        is accessed by FTP and such file is downloaded. This does not
        require additional libraries. However in order
        to produce FASTA of FASTQ files one would need to use SRA-Toolkit.
        Thus, it is assumed that this library is already installed or it
        will be installed in the near future. One can immediately specify
        the download type to fasta or fastq.

        Following  ``**kwargs`` can be passed:
            * filetype - str
                can be sra, fasta, or fastq - for fasta or fastq SRA-Toolkit
                need to be installed
            * aspera - bool
                use Aspera to download samples, defaults to False
            * keep_sra - bool
                keep SRA files after download. Removes SRA file only if the
                selected file type is different than "sra", defaults to False
            * fastq_dump_options - dict
                pass options to fastq-dump (if used, the options has to be in
                long form eg. --split-files), defaults to::
                    {
                        'split-files': None,
                        'readids': None,
                        'read-filter': 'pass',
                        'dumpbase': None,
                        'gzip': None
                    }

        Args:
            email (:obj:`str`): an email (any) - Required by NCBI for access
            directory (:obj:`str`, optional): The directory to which download
                the data. Defaults to "./".
            **kwargs: Arbitrary keyword arguments, see description

        Returns:
            :obj:`list` of :obj:`str`: List of downloaded files.

        Raises:
            :obj:`TypeError`: Type to download unknown
            :obj:`NoSRARelationException`: No SRAToolkit
            :obj:`Exception`: Wrong e-mail
            :obj:`HTTPError`: Cannot access or connect to DB
        """
        from Bio import Entrez

        # Unpack arguments
        filetype = kwargs.get('filetype', 'fasta')
        aspera = kwargs.get('aspera', False)
        keep_sra = kwargs.get('keep_sra', False)

        fq_options = {
            'split-files': None,
            'readids': None,
            'read-filter': 'pass',
            'dumpbase': None,
            'gzip': None
        }
        fastq_dump_options = kwargs.get('fastq_dump_options', fq_options)

        # Check download filetype
        filetype = filetype.lower()
        if filetype not in ["sra", "fastq", "fasta"]:
            raise TypeError("Unknown type to downlod: %s. Use sra, fastq or "
                            "fasta." % filetype)

        # Setup the query
        ftpaddres = ("ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/"
                     "ByExp/sra/SRX/{range_subdir}/{record_dir}/{file_dir}/"
                     "{file_dir}.sra")
        queries = []
        try:
            for sra in self.relations['SRA']:
                query = sra.split("=")[-1]
                assert 'SRX' in query, "Sample looks like it is not SRA: %s" % query
                logger.info("Query: %s" % query)
                queries.append(query)
        except KeyError:
            raise NoSRARelationException(
                'No relation called SRA for %s' % self.get_accession())

        # check if the e-mail is more or less not a total crap
        Entrez.email = email
        if not (Entrez.email is not None and '@' in email and
                        email != '' and '.' in email):
            raise Exception('You have to provide valid e-mail')

        downloaded_paths = list()
        for query in queries:
            # retrieve IDs for given SRX
            searchdata = Entrez.esearch(db='sra', term=query, usehistory='y',
                                        retmode='json')
            answer = json.loads(searchdata.read())
            ids = answer["esearchresult"]["idlist"]
            assert len(ids) == 1, "There should be one and only one ID per SRX"

            # using ID fetch the info
            number_of_trials = 10
            wait_time = 30
            for trial in range(number_of_trials):
                try:
                    results = Entrez.efetch(db="sra", id=ids[0],
                                            rettype="runinfo",
                                            retmode="text").read()
                    break
                except HTTPError as httperr:
                    if "502" in str(httperr):
                        logger.error(("Error: %s, trial %i out of %i, waiting "
                                      "for %i seconds.") % (
                                         str(httperr),
                                         trial,
                                         number_of_trials,
                                         wait_time))
                        time.sleep(wait_time)
                    else:
                        raise httperr
            try:
                df = DataFrame(
                    [i.split(',') for i in results.split('\n') if i != ''][1:],
                    columns=
                    [i.split(',') for i in results.split('\n') if i != ''][0])
            except IndexError:
                logger.error(("SRA is empty (ID: %s, query: %s). "
                              "Check if it is publicly available.") %
                             (ids[0], query))
                continue

            # check it first
            try:
                df['download_path']
            except KeyError as e:
                logger.error('KeyError: ' + str(e) + '\n')
                logger.error(str(results) + '\n')

            # make the directory
            if platform.system() == "Windows":
                name_regex = r'[\s\*\?\(\),\.\:\%\|\"\<\>]'
            else:
                name_regex = r'[\s\*\?\(\),\.;]'
            directory_path = os.path.abspath(
                os.path.join(directory, "%s_%s_%s" % (
                    'Supp',
                    self.get_accession(),
                    # the directory name cannot contain many of the signs
                    re.sub(name_regex, '_', self.metadata['title'][0]))))

            utils.mkdir_p(os.path.abspath(directory_path))
            for path in df['download_path']:
                sra_run = path.split("/")[-1]
                logger.info("Analysing %s" % sra_run)
                url = ftpaddres.format(range_subdir=query[:6],
                                       record_dir=query,
                                       file_dir=sra_run)
                logger.debug("URL: %s", url)
                filepath = os.path.abspath(
                    os.path.join(directory_path, "%s.sra" % sra_run))
                utils.download_from_url(url, filepath, aspera=aspera)

                if filetype in ["fasta", "fastq"]:
                    if utils.which('fastq-dump') is None:
                        raise NoSRAToolkitException(
                            "fastq-dump command not found")
                    ftype = ""
                    if filetype == "fasta":
                        ftype = " --fasta "
                    cmd = "fastq-dump"
                    for fqoption, fqvalue in iteritems(fastq_dump_options):
                        if fqvalue:
                            cmd += (" --%s %s" % (fqoption, fqvalue))
                        else:
                            cmd += (" --%s" % fqoption)
                    cmd += " %s --outdir %s %s"
                    cmd = cmd % (ftype, directory_path, filepath)
                    logger.debug(cmd)
                    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE,
                                               shell=True)
                    logger.info("Converting to %s/%s*.%s.gz\n" % (
                        directory_path, sra_run, filetype))
                    pout, perr = process.communicate()
                    downloaded_paths += glob.glob(os.path.join(
                        directory_path,
                        "%s*.%s.gz" % (sra_run, filetype)
                    ))
                    if not keep_sra and filetype != 'sra':
                        # Delete sra file
                        os.unlink(filepath)
                    else:
                        downloaded_paths.append(filepath)

        return downloaded_paths


class GPL(SimpleGEO):
    """Class that represents platform from GEO database"""

    geotype = "PLATFORM"

    def __init__(self, name, metadata, table=None, columns=None, gses=None,
                 gsms=None, database=None):
        """Initialize GPL.

        Args:
            name (:obj:`str`): Name of the object
            metadata (:obj:`dict`): Metadata information
            table (:obj:`pandas.DataFrame`, optional): Table with actual GPL data
            columns (:obj:`pandas.DataFrame`, optional): Table with description
                of the columns. Defaults to None.
            gses (:obj:`dict` of :obj:`GEOparse.GSE`, optional): A dictionary of
                GSE objects. Defaults to None.
            gsms (:obj:`dict` of :obj:`GEOparse.GSM`, optional): A dictionary of
                GSM objects. Defaults to None.
            database (:obj:`GEOparse.GEODatabase`, optional): A database object
                from SOFT file associated with GPL. Defaults to None.
        """
        gses = {} if gses is None else gses
        if not isinstance(gses, dict):
            raise ValueError(
                "GSEs should be a dictionary not a %s" % str(type(gses)))
        gsms = {} if gsms is None else gsms
        if not isinstance(gsms, dict):
            raise ValueError(
                "GSMs should be a dictionary not a %s" % str(type(gsms)))

        for gsm_name, gsm in iteritems(gsms):
            assert isinstance(gsm, GSM), "All GSMs should be of type GSM"
        for gse_name, gse in iteritems(gses):
            assert isinstance(gse, GSE), "All GSEs should be of type GSE"
        if database is not None:
            if not isinstance(database, GEODatabase):
                raise ValueError(
                    "Database should be a GEODatabase not a %s" % str(
                        type(database)))

        table = DataFrame() if table is None else table
        columns = DataFrame() if columns is None else columns
        SimpleGEO.__init__(self, name=name, metadata=metadata, table=table,
                           columns=columns)

        self.gses = gses
        self.gsms = gsms
        self.database = database


class GDSSubset(BaseGEO):
    """Class that represents a subset from GEO GDS object."""

    geotype = "SUBSET"

    def _get_object_as_soft(self):
        """Get the object as SOFT formatted string."""
        soft = ["^%s = %s" % (self.geotype, self.name),
                self._get_metadata_as_string()]
        return "\n".join(soft)


class GEODatabase(BaseGEO):
    """Class that represents a subset from GEO GDS object."""

    geotype = "DATABASE"

    def _get_object_as_soft(self):
        """Return object as SOFT formatted string."""
        soft = ["^%s = %s" % (self.geotype, self.name),
                self._get_metadata_as_string()]
        return "\n".join(soft)


class GDS(SimpleGEO):
    """Class that represents a dataset from GEO database"""

    geotype = "DATASET"

    def __init__(self, name, metadata, table, columns, subsets, database=None):
        """Initialize GDS

        Args:
            name (:obj:`str`): Name of the object.
            metadata (:obj:`dict`): Metadata information.
            table (:obj:`pandas.DataFrame`): Table with the data from SOFT file.
            columns (:obj:`pandas.DataFrame`): description of the columns,
                number of columns, order, and names represented as index in
                this DataFrame has to be the same as table.columns.
            subsets (:obj:`dict` of :obj:`GEOparse.GDSSubset`): GDSSubset from
                GDS soft file.
            database (:obj:`GEOparse.Database`, optional): Database from SOFT
                file. Defaults to None.
        """
        if not isinstance(subsets, dict):
            raise ValueError(
                "Subsets should be a dictionary not a %s" % str(type(subsets)))
        if database is not None:
            if not isinstance(database, GEODatabase):
                raise ValueError(
                    "Database should be a GEODatabase not a %s" % str(
                        type(database)))

        SimpleGEO.__init__(self, name=name, metadata=metadata, table=table,
                           columns=columns)
        # effectively deletes the columns with ID_REF
        self.columns = self.columns.dropna()
        self.subsets = subsets
        self.database = database

        for subset_name, subset in iteritems(subsets):
            message = "All subsets should be of type GDSSubset"
            assert isinstance(subset, GDSSubset), message

    def _get_object_as_soft(self):
        """Return object as SOFT formatted string."""
        soft = []
        if self.database is not None:
            soft.append(self.database._get_object_as_soft())
        soft += ["^%s = %s" % (self.geotype, self.name),
                 self._get_metadata_as_string()]
        for subset in self.subsets.values():
            soft.append(subset._get_object_as_soft())
        soft += ["^%s = %s" % (self.geotype, self.name),
                 self._get_columns_as_string(),
                 self._get_table_as_string()]
        return "\n".join(soft)


class GSE(BaseGEO):
    """Class representing GEO series"""

    geotype = "SERIES"

    def __init__(self, name, metadata, gpls=None, gsms=None, database=None):
        """Initialize GSE.

        Args:
            name (:obj:`str`): Name of the object.
            metadata (:obj:`dict`): Metadata information.
            gpls (:obj:`dict` of :obj:`GEOparse.GPL`, optional): A dictionary of
                GSE objects. Defaults to None.
            gsms (:obj:`dict` of :obj:`GEOparse.GSM`, optional): A dictionary of
                GSM objects. Defaults to None.
            database (:obj:`GEOparse.Database`, optional): Database from SOFT
                file. Defaults to None.
        """
        gpls = {} if gpls is None else gpls
        if not isinstance(gpls, dict):
            raise ValueError(
                "GPLs should be a dictionary not a %s" % str(type(gpls)))
        gsms = {} if gsms is None else gsms
        if not isinstance(gsms, dict):
            raise ValueError(
                "GSMs should be a dictionary not a %s" % str(type(gsms)))

        for gsm_name, gsm in iteritems(gsms):
            assert isinstance(gsm, GSM), "All GSMs should be of type GSM"
        for gpl_name, gpl in iteritems(gpls):
            assert isinstance(gpl, GPL), "All GPLs should be of type GPL"
        if database is not None:
            if not isinstance(database, GEODatabase):
                raise ValueError(
                    "Database should be a GEODatabase not a %s" % str(
                        type(database)))
        BaseGEO.__init__(self, name=name, metadata=metadata)

        self.gpls = gpls
        self.gsms = gsms
        self.database = database
        self._phenotype_data = None

    @property
    def phenotype_data(self):
        """Get the phenotype data for each of the sample."""
        if self._phenotype_data is None:
            pheno_data = {}
            for gsm_name, gsm in iteritems(self.gsms):
                tmp = {}
                for key, value in iteritems(gsm.metadata):
                    if len(value) == 0:
                        tmp[key] = np.nan
                    elif key.startswith("characteristics_"):
                        for i, char in enumerate(value):
                            char = re.split(":\s+", char)
                            char_type, char_value = [char[0],
                                                     ": ".join(char[1:])]
                            tmp[key + "." + str(
                                i) + "." + char_type] = char_value
                    else:
                        tmp[key] = ",".join(value)
                pheno_data[gsm_name] = tmp
            self._phenotype_data = DataFrame(pheno_data).T
        return self._phenotype_data

    def merge_and_average(self, platform, expression_column, group_by_column,
                          force=False, merge_on_column=None, gsm_on=None,
                          gpl_on=None):
        """Merge and average GSE samples.

        For given platform prepare the DataFrame with all the samples present in
        the GSE annotated with given column from platform and averaged over
        the column.

        Args:
            platform (:obj:`str` or :obj:`GEOparse.GPL`): GPL platform to use.
            expression_column (:obj:`str`): Column name in which "expressions"
                are represented
            group_by_column (:obj:`str`): The data will be grouped and averaged
                over this column and only this column will be kept
            force (:obj:`bool`): If the name of the GPL does not match the
                platform name in GSM proceed anyway
            merge_on_column (:obj:`str`): Column to merge the data on - should
                be present in both GSM and GPL
            gsm_on (:obj:`str`): In the case columns to merge are different in
                GSM and GPL use this column in GSM
            gpl_on (:obj:`str`): In the case columns to merge are different in
                GSM and GPL use this column in GPL

        Returns:
            :obj:`pandas.DataFrame`: Merged and averaged table of results.

        """
        if isinstance(platform, str):
            gpl = self.gpls[platform]
        elif isinstance(platform, GPL):
            gpl = platform
        else:
            raise ValueError("Platform has to be of type GPL or string with "
                             "key for platform in GSE")

        data = []
        for gsm in self.gsms.values():
            if gpl.name == gsm.metadata['platform_id'][0]:
                data.append(gsm.annotate_and_average(
                    gpl=gpl,
                    merge_on_column=merge_on_column,
                    expression_column=expression_column,
                    group_by_column=group_by_column,
                    force=force,
                    gpl_on=gpl_on,
                    gsm_on=gsm_on))
        if len(data) == 0:
            logger.warning("No samples for the platform were found\n")
            return None
        elif len(data) == 1:
            return data[0]
        else:
            return data[0].join(data[1:])

    def pivot_samples(self, values, index="ID_REF"):
        """Pivot samples by specified column.

        Construct a table in which columns (names) are the samples, index
        is a specified column eg. ID_REF and values in the columns are of one
        specified type.

        Args:
            values (:obj:`str`): Column name present in all GSMs.
            index (:obj:`str`, optional): Column name that will become an index in
                pivoted table. Defaults to "ID_REF".

        Returns:
            :obj:`pandas.DataFrame`: Pivoted data

        """
        data = []
        for gsm in self.gsms.values():
            tmp_data = gsm.table.copy()
            tmp_data["name"] = gsm.name
            data.append(tmp_data)
        ndf = concat(data).pivot(index=index, values=values, columns="name")
        return ndf

    def pivot_and_annotate(self, values, gpl, annotation_column, gpl_on="ID",
                           gsm_on="ID_REF"):
        """Annotate GSM with provided GPL.

        Args:
            values (:obj:`str`): Column to use as values eg. "VALUES"
            gpl (:obj:`pandas.DataFrame` or :obj:`GEOparse.GPL`): A Platform or
                DataFrame to annotate with.
            annotation_column (:obj:`str`): Column in table for annotation.
            gpl_on (:obj:`str`, optional): Use this column in GPL to merge.
                Defaults to "ID".
            gsm_on (:obj:`str`, optional): Use this column in GSM to merge.
                Defaults to "ID_REF".

        Returns:
            pandas.DataFrame: Pivoted and annotated table of results

        """
        if isinstance(gpl, GPL):
            annotation_table = gpl.table
        elif isinstance(gpl, DataFrame):
            annotation_table = gpl
        else:
            raise TypeError("gpl should be a GPL object or a pandas.DataFrame")
        pivoted_samples = self.pivot_samples(values=values, index=gsm_on)
        ndf = pivoted_samples.reset_index().merge(
            annotation_table[[gpl_on, annotation_column]],
            left_on=gsm_on,
            right_on=gpl_on).set_index(gsm_on)
        del ndf[gpl_on]
        ndf.columns.name = 'name'
        return ndf

    def download_supplementary_files(self, directory='series',
                                     download_sra=True, sra_filetype='fasta',
                                     email=None, sra_kwargs=None):
        """Download supplementary data.

        Args:
            directory (:obj:`str`, optional): Directory to download the data
                (in this directory function will create new directory with the
                files), by default this will be named with the series
                name + _Supp.
            download_sra (:obj:`bool`, optional): Indicates whether to download
                SRA raw data too. Defaults to True.
            sra_filetype (:obj:`str`, optional): Indicates what file type to
                download if we specified SRA, can be sra, fasta or fastq.
                Defaults to "fasta".
            email (:obj:`str`, optional): E-mail that will be provided to the
                Entrez. Defaults to None.
            sra_kwargs (:obj:`dict`, optional): Kwargs passed to the
                GSM.download_SRA method. Defaults to None.

        Returns:
            :obj:`dict`: Downloaded paths for each of the GSM
        """
        if directory == 'series':
            dirpath = os.path.abspath(self.get_accession() + "_Supp")
            utils.mkdir_p(dirpath)
        else:
            dirpath = os.path.abspath(directory)
            utils.mkdir_p(dirpath)
        downloaded_paths = dict()
        for gsmname, gsm in iteritems(self.gsms):
            paths = gsm.download_supplementary_files(email=email,
                                                     download_sra=download_sra,
                                                     directory=dirpath,
                                                     sra_kwargs=sra_kwargs)
            downloaded_paths[gsmname] = paths

        return downloaded_paths

    def download_SRA(self, email, directory='series', filterby=None, **kwargs):
        """Download SRA files for each GSM in series.

        Args:
            email (:obj:`str`): E-mail that will be provided to the Entrez.
            directory (:obj:`str`, optional): Directory to save the data
                (defaults to the 'series' which saves the data to the directory
                with the name of the series + '_SRA' ending).
                Defaults to "series".
            filterby (:obj:`str`, optional): Filter GSM objects, argument is a
                function that operates on GSM object  and return bool
                eg. lambda x: "brain" not in x.name. Defaults to None.
            **kwargs: Any arbitrary argument passed to GSM.download_SRA
                method. See the documentation for more details.
        """
        if directory == 'series':
            dirpath = os.path.abspath(self.get_accession() + "_SRA")
            utils.mkdir_p(dirpath)
        else:
            dirpath = os.path.abspath(directory)
            utils.mkdir_p(dirpath)
        if filterby is not None:
            gsms_to_use = [gsm for gsm in self.gsms.values() if filterby(gsm)]
        else:
            gsms_to_use = self.gsms.values()

        downloaded_paths = dict()
        for gsm in gsms_to_use:
            logger.info(
                "Downloading SRA files for %s series\n" % (gsm.name))
            paths = gsm.download_SRA(email=email, directory=dirpath, **kwargs)
            downloaded_paths[gsm.name] = paths

        return downloaded_paths

    def _get_object_as_soft(self):
        """Get object as SOFT formatted string."""
        soft = []
        if self.database is not None:
            soft.append(self.database._get_object_as_soft())
        soft += ["^%s = %s" % (self.geotype, self.name),
                 self._get_metadata_as_string()]
        for gsm in itervalues(self.gsms):
            soft.append(gsm._get_object_as_soft())
        for gpl in itervalues(self.gpls):
            soft.append(gpl._get_object_as_soft())

        return "\n".join(soft)

    def __str__(self):
        return str("<%s: %s - %i SAMPLES, %i PLATFORM(s)>" % (
            self.geotype, self.name, len(self.gsms), len(self.gpls)))

    def __repr__(self):
        return str("<%s: %s - %i SAMPLES, %i PLATFORM(s)>" % (
            self.geotype, self.name, len(self.gsms), len(self.gpls)))
