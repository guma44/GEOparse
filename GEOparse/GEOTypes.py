"""
Classes that represent different GEO entities
"""

import os
import re
import abc
import sys
import gzip
import json
import time
import subprocess
import numpy as np
from . import utils
from sys import stderr, stdout
from pandas import DataFrame, concat
try:
    from urllib.error import HTTPError
except ImportError:
    from urllib2 import HTTPError
from six import iteritems, itervalues



class DataIncompatibilityException(Exception): pass
class NoMetadataException(Exception): pass
class NoSRARelationException(Exception): pass
class NoSRAToolkitException(Exception): pass


class BaseGEO(object):

    __metaclass__ = abc.ABCMeta

    geotype = None

    def __init__(self, name, metadata):
        """base GEO object

        :param name: str -- name of the object
        :param metadata: dict -- metadata information
        """

        if not isinstance(metadata, dict):
            raise ValueError("Metadata should be a dictionary not a %s" % str(type(metadata)))

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

        :param metaname: str -- name of the attribute
        :returns: list or str -- if there is more than one attribute it returns a list
        otherwise it returns a string

        """
        metadata_value = self.metadata.get(metaname, None)
        if metadata_value is None:
            raise NoMetadataException("No metadata attribute named %s" % metaname)
        if not isinstance(metadata_value, list):
            raise TypeError("Metadata is not a list and it should be.")

        if len(metadata_value) > 1:
            return metadata_value
        else:
            return metadata_value[0]

    def get_accession(self):
        """Return accession ID of the sample
        :returns: str

        """
        return self.get_metadata_attribute("geo_accession")

    def get_type(self):
        """Return type of a subset
        """
        try:
            return self.get_metadata_attribute("type")
        except NoMetadataException:
            return None

    def _get_metadata_as_string(self):
        """Returns metadata as SOFT formated string
        """
        metalist = []
        for metaname, meta in iteritems(self.metadata):
            assert isinstance(meta, list), "Single value in metadata dictionary should be a list!"
            for data in meta:
                if data:
                    metalist.append("!%s_%s = %s" % (self.geotype.capitalize(),
                                                     metaname, data))
        return "\n".join(metalist)

    def show_metadata(self):
        """
         Show metadat in SOFT format
        """
        print(self._get_metadata_as_string())

    def to_soft(self, path_or_handle, as_gzip=False):
        """Save the object in a SOFT format.

        :param path_or_handle: path or handle to output file
        :param gzip: bool -- save as gzip
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
        """
         Return object as SOFT formated string.
        """
        raise NotImplementedError("Method not implemented")

    def __str__(self):
        return str("<%s: %s>" % (self.geotype, self.name))

    def __repr__(self):
        return str("<%s: %s>" % (self.geotype, self.name))


class SimpleGEO(BaseGEO):

    __metaclass__ = abc.ABCMeta

    def __init__(self, name, metadata, table, columns):
        """base GEO object

        :param name: str -- name of the object
        :param table: pandas.DataFrame -- table with the data from SOFT file
        :param metadata: dict -- metadata information
        :param columns: pandas.DataFrame -- description of the columns, number of columns, order, and names
        represented as index in this DataFrame has to be the same as table.columns.

        """
        if not isinstance(table, DataFrame):
            raise ValueError("Table data should be an instance of pandas.DataFrame not %s" % str(type(table)))
        if not isinstance(columns, DataFrame):
            raise ValueError("Columns description should be an instance of pandas.DataFrame not %s" % str(type(columns)))

        BaseGEO.__init__(self, name=name, metadata=metadata)

        self.table = table
        self.columns = columns
        if self.columns.index.tolist() != self.table.columns.tolist():
            if sorted(self.columns.index.tolist()) == sorted(self.table.columns.tolist()):
                stderr.write("Data columns in %s %s are not in order. Reordering.\n" % (self.geotype, self.name))
                self.columns = self.columns.ix[self.table.columns]
            else:
                rows_in_columns = ", ".join(self.columns.index.tolist())
                columns_in_table = ", ".join(self.table.columns.tolist())
                raise DataIncompatibilityException("\nData columns do not match columns description index in %s\n" % (self.name) +
                                                   "Columns in table are: %s\n" % columns_in_table +
                                                   "Index in columns are: %s\n" % rows_in_columns
                                                   )
        if self.columns.columns[0] != 'description':
            raise ValueError("Columns table must contain a column named 'description'. Here columns are: %s" % ", ".join(map(str, self.columns.columns)))

    def head(self):
        """Return short description of the object
        """
        stdout.write("%s %s" % (self.geotype, self.name) + "\n")
        stdout.write(" - Metadata:" + "\n")
        stdout.write("\n".join(self._get_metadata_as_string().split("\n")[:5]) + "\n")
        stdout.write("\n")
        stdout.write(" - Columns:" + "\n")
        stdout.write(self.columns.to_string() + "\n")
        stdout.write("\n")
        stdout.write(" - Table:" + "\n")
        stdout.write("\t".join(["Index"] + self.table.columns.tolist()) + "\n")
        stdout.write(self.table.head().to_string(header=None) + "\n")
        stdout.write(" "*40 + "..." + " "*40 + "\n")
        stdout.write(" "*40 + "..." + " "*40 + "\n")
        stdout.write(" "*40 + "..." + " "*40 + "\n")
        stdout.write(self.table.tail().to_string(header=None) + "\n")

    def show_columns(self):
        """Show columns in SOFT format
        """
        print(self.columns)

    def show_table(self, number_of_lines=5):
        """
        Show few lines of the table the table as pandas.DataFrame

        :param number_of_lines: int -- number of lines to show
        """
        print(self.table.head(number_of_lines))

    def _get_object_as_soft(self):
        """
         Return object as SOFT formated string.
        """
        soft = ["^%s = %s" % (self.geotype, self.name),
                self._get_metadata_as_string(),
                self._get_columns_as_string(),
                self._get_table_as_string()]
        return "\n".join(soft)

    def _get_table_as_string(self):
        """Returns table as SOFT formated string
        """
        tablelist = []
        tablelist.append("!%s_table_begin" % self.geotype.lower())
        tablelist.append("\t".join(self.table.columns))
        for idx, row in self.table.iterrows():
            tablelist.append("\t".join(map(str, row)))
        tablelist.append("!%s_table_end" % self.geotype.lower())
        return "\n".join(tablelist)

    def _get_columns_as_string(self):
        """Returns columns  as SOFT formated string
        """
        columnslist = []
        for rowidx, row in self.columns.iterrows():
            columnslist.append("#%s = %s" % (rowidx, row.description))
        return "\n".join(columnslist)


class GSM(SimpleGEO):

    """Class that represents sample from GEO database"""

    geotype = 'SAMPLE'

    def annotate(self, gpl, annotation_column, gpl_on="ID", gsm_on="ID_REF", in_place=False):
        """Annotate GSM with provided GPL

        :param gpl: pandas.DataFrame or GPL -- a Platform or DataFrame to annotate with
        :param annotation_column: str -- column in table for annotation
        :param gsm_on: str -- use this column in GSM to merge, defaults to ID_REF
        :param gpl_on: str -- use this column in GPL to merge, defaults to ID
        :param in_place: bool -- if True substitute table in GSM by new annotated table, defaults to False
        :returns: DataFrame or None if in_place=True

        """
        if isinstance(gpl, GPL):
            annotation_table = gpl.table
        elif isinstance(gpl, DataFrame):
            annotation_table = gpl
        else:
            raise TypeError("gpl should be a GPL object or a pandas.DataFrame")

        # annotate by merging
        annotated = self.table.merge(annotation_table[[gpl_on, annotation_column]], left_on=gsm_on, right_on=gpl_on)
        del annotated[gpl_on]
        if in_place:
            self.table = annotated
            return None
        else:
            return annotated

    def annotate_and_average(self, gpl, expression_column, group_by_column, rename=True,
                             force=False, merge_on_column=None, gsm_on=None, gpl_on=None):
        """Annotate GSM table with provided GPL

        :param gpl: GPL object -- platform for annotations
        :param expression_column: str -- column name in which "expressions" are represented
        :param group_by_column: str -- the data will be grouped and averaged over this column and only this column will be kept
        :param rename: bool -- rename output column to the self.name
        :param force: bool -- if the name of the GPL does not match the platform name in GSM proceed anyway
        :param merge_on_column: str -- column to merge the data on - should be present in both GSM and GPL
        :param gsm_on: str -- in the case columns to merge are different in GSM and GPL use this column in GSM
        :param gpl_on: str -- in the case columns to merge are different in GSM and GPL use this column in GPL
        :returns: pandas.DataFrame

        """
        if gpl.name != self.metadata['platform_id'][0] and not force:
            raise KeyError("Platforms from GSM (%s) and from GPL (%s)" % (gpl.name, self.metadata['platform_id']) +
                           " are incompatible. Use force=True to use this GPL.")
        if merge_on_column is None and gpl_on is None and gsm_on is None:
            raise Exception("You have to provide one of the two: merge_on_column or gpl_on and gsm_on parameters")
        if merge_on_column:
            stderr.write("merge_on_column is not None. Using this option.\n")
            tmp_data = self.table.merge(gpl.table, on=merge_on_column, how='outer')
            tmp_data = tmp_data.groupby(group_by_column).mean()[[expression_column]]
        else:
            if gpl_on is None or gsm_on is None:
                raise Exception("Please provide both gpl_on and gsm_on or provide merge_on_column only")
            tmp_data = self.table.merge(gpl.table, left_on=gsm_on, right_on=gpl_on, how='outer')
            tmp_data = tmp_data.groupby(group_by_column).mean()[[expression_column]]
        if rename:
            tmp_data.columns = [self.name]
        return tmp_data

    def download_supplementary_files(self, directory="./", download_sra=True, sra_filetype='fasta', email=None):
        """Download all supplementary data available for the sample

        :param directory: directory to download the data (in this directory function will create
                          new directory with the files)
        :param download_sra: bool - indicates whether to download SRA raw data too, defaults to True
        :param sra_filetype: indicates what file type to download if we specified SRA, can be sra, fasta or fastq
        :param email: e-mail that will be provided to the Entrez, defaults to None
        """
        directory_path = os.path.abspath(os.path.join(directory, "%s_%s_%s" % ('Supp',
                                                                               self.get_accession(),
                                                                               re.sub(r'[\s\*\?\(\),\.;]', '_', self.metadata['title'][0]) # the directory name cannot contain many of the signs
                                                                               )))
        utils.mkdir_p(os.path.abspath(directory_path))
        for metakey, metavalue in iteritems(self.metadata):
            if 'supplementary_file' in metakey:
                assert len(metavalue) == 1 and metavalue != ''
                # stderr.write("Downloading %s\n" % metavalue)
                if 'sra' in metavalue[0] and download_sra:
                    self.download_SRA(email, filetype=sra_filetype, directory=directory)
                else:
                    download_path = os.path.abspath(os.path.join(directory, os.path.join(directory_path, metavalue[0].split("/")[-1])))
                    utils.download_from_url(metavalue[0], download_path)



    def download_SRA(self, email, metadata_key='auto', directory='./', filetype='sra', aspera=False, keep_sra=False):
        """Download RAW data as SRA file to the sample directory created ad hoc
        or the directory specified by the parameter. The sample has to come from
        sequencing eg. mRNA-seq, CLIP etc.

        An important parameter is a download_type. By default an SRA is accessed by FTP and
        such file is downloaded. This does not require additional libraries. However in order
        to produce FASTA of FASTQ files one would need to use SRA-Toolkit. Thus, it is assumed
        that this library is already installed or it will be installed in the near future. One
        can immediately specify the download type to fasta or fastq.

        :param email: an email (any) - required by NCBI for access
        :param directory: The directory to which download the data. By default current directory is used
        :param filetype: can be sra, fasta, or fastq - for fasta or fastq SRA-Toolkit need to be installed
        :param aspera: bool - use Aspera to download samples, defaults to False
        :param keep_sra: bool - keep SRA files after download, defaults to False

        """
        from Bio import Entrez
        # Check download filetype
        filetype = filetype.lower()
        if filetype not in ["sra", "fastq", "fasta"]:
            raise Exception("Unknown type to downlod: %s. Use sra, fastq or fasta." % filetype)

        # Setup the query
        ftpaddres = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/{range_subdir}/{record_dir}/{file_dir}/{file_dir}.sra"
        queries = []
        try:
            for sra in self.relations['SRA']:
                query = sra.split("=")[-1]
                assert 'SRX' in query, "Sample looks like it is not SRA: %s" % query
                print("Query: %s" % query)
                queries.append(query)
        except KeyError:
            raise NoSRARelationException('No relation called SRA for %s' % self.get_accession())

        # check if the e-mail is more or less not a total crap
        Entrez.email = email
        if not (Entrez.email is not None and '@' in email and email != '' and '.' in email):
            raise Exception('You have to provide valid e-mail')

        for query in queries:
            # retrieve IDs for given SRX
            searchdata = Entrez.esearch(db='sra', term=query, usehistory='y', retmode='json')
            answer = json.loads(searchdata.read())
            ids = answer["esearchresult"]["idlist"]
            assert len(ids) == 1, "There should be one and only one ID per SRX"

            # using ID fetch the info
            number_of_trials = 10
            wait_time = 30
            for tiral in range(number_of_trials):
                try:
                    results = Entrez.efetch(db="sra", id=ids[0], rettype="runinfo", retmode="text").read()
                    break
                except HTTPError as httperr:
                    if "502" in str(httperr):
                        sys.stderr.write("Error: %s, trial %i out of %i, waiting for %i seconds." % (str(httperr),
                                                                                                     trial,
                                                                                                     number_of_trials,
                                                                                                     wait_time))
                        time.sleep(wait_time)
                    else:
                        raise httperr
            df = DataFrame([i.split(',') for i in results.split('\n') if i != ''][1:], columns = [i.split(',') for i in results.split('\n') if i != ''][0])

            # check it first
            try:
                df['download_path']
            except KeyError as e:
                stderr.write('KeyError: ' + str(e) + '\n')
                stderr.write(str(results) + '\n')

            # make the directory
            if platform.system() == "Windows":
                name_regex = r'[\s\*\?\(\),\.\:\%\|\"\<\>]'
            else:
                name_regex = r'[\s\*\?\(\),\.;]'
            directory_path = os.path.abspath(os.path.join(directory, "%s_%s_%s" % ('Supp',
                                                                                   self.get_accession(),
                                                                                   re.sub(name_regex, '_', self.metadata['title'][0]) # the directory name cannot contain many of the signs
                                                                                   )))
            utils.mkdir_p(os.path.abspath(directory_path))

            for path in df['download_path']:
                sra_run = path.split("/")[-1]
                print("Analysing %s" % sra_run)
                url = ftpaddres.format(range_subdir=query[:6],
                                           record_dir=query,
                                           file_dir=sra_run)
                filepath = os.path.abspath(os.path.join(directory_path, "%s.sra" % sra_run))
                utils.download_from_url(url, filepath, aspera=aspera)

                if filetype in ["fasta", "fastq"]:
                    ftype = ""
                    if filetype == "fasta":
                        ftype = " --fasta "
                    cmd = "fastq-dump --split-files --gzip %s --outdir %s %s"
                    cmd = cmd % (ftype, directory_path, filepath)

                    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    stderr.write("Converting to %s/%s_*.%s.gz\n" % (
                        directory_path, sra_run, filetype))
                    pout, perr = process.communicate()
                    if "command not found" in perr:
                        raise NoSRAToolkitException("fastq-dump command not found")
                    else:
                        print(pout)
                        if not keep_sra:
                            # Delete sra file
                            os.unlink(filepath)

class GPL(SimpleGEO):

    """Class that represents platform from GEO database"""

    geotype = "PLATFORM"

    def __init__(self, name, metadata, table=None, columns=None, gses=None, gsms=None,  database=None):
        """Initialize GPL

        :param name: str -- name of the object
        :param metadata: dict -- metadata information
        :param gses: list -- list of GSE objects
        :param gsms: list -- list of GSM objects
        :param database: GEODatabase -- Database from SOFT file
        """

        gses = {} if gses is None else gses
        if not isinstance(gses, dict):
            raise ValueError("GSEs should be a dictionary not a %s" % str(type(gses)))
        gsms = {} if gsms is None else gsms
        if not isinstance(gsms, dict):
            raise ValueError("GSMs should be a dictionary not a %s" % str(type(gsms)))

        for gsm_name, gsm in iteritems(gsms):
            assert isinstance(gsm, GSM), "All GSMs should be of type GSM"
        for gse_name, gse in iteritems(gses):
            assert isinstance(gse, GSE), "All GSEs should be of type GSE"
        if database is not None:
            if not isinstance(database, GEODatabase):
                raise ValueError("Database should be a GEODatabase not a %s" % str(type(database)))

        table = DataFrame() if table is None else table
        columns = DataFrame() if columns is None else columns
        SimpleGEO.__init__(self, name=name, metadata=metadata, table=table, columns=columns)

        self.gses = gses
        self.gsms = gsms
        self.database = database

class GDSSubset(BaseGEO):

    """Class that represents a subset from GEO GDS object"""

    geotype = "SUBSET"

    def _get_object_as_soft(self):
        """
         Return object as SOFT formated string.
        """
        soft = ["^%s = %s" % (self.geotype, self.name),
                self._get_metadata_as_string()]
        return "\n".join(soft)


class GEODatabase(BaseGEO):

    """Class that represents a subset from GEO GDS object"""

    geotype = "DATABASE"

    def _get_object_as_soft(self):
        """
         Return object as SOFT formated string.
        """
        soft = ["^%s = %s" % (self.geotype, self.name),
                self._get_metadata_as_string()]
        return "\n".join(soft)


class GDS(SimpleGEO):

    """Class that represents a dataset from GEO database"""

    geotype = "DATASET"

    def __init__(self, name, metadata, table, columns, subsets, database=None):
        """Initialize GDS

        :param name: str -- name of the object
        :param table: pandas.DataFrame -- table with the data from SOFT file
        :param metadata: dict -- metadata information
        :param columns: pandas.DataFrame -- description of the columns, number of columns, order, and names
        represented as index in this DataFrame has to be the same as table.columns.
        :param subsets: dict -- dictionary of GDSSubset from GDS soft file
        :param database: GEODatabase -- Database from SOFT file
        """
        if not isinstance(subsets, dict):
            raise ValueError("Subsets should be a dictionary not a %s" % str(type(subsets)))
        if database is not None:
            if not isinstance(database, GEODatabase):
                raise ValueError("Database should be a GEODatabase not a %s" % str(type(database)))


        SimpleGEO.__init__(self, name=name, metadata=metadata, table=table, columns=columns)
        self.columns = self.columns.dropna() # effectively deletes the columns with ID_REF
        self.subsets = subsets
        self.database = database

        for subset_name, subset in iteritems(subsets):
            assert isinstance(subset, GDSSubset), "All subsets should be of type GDSSubset"

    def _get_object_as_soft(self):
        """
         Return object as SOFT formated string.
        """
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
        """Initialize GSE

        :param name: str -- name of the object
        :param metadata: dict -- metadata information
        :param gpls: list -- list of GPL objects
        :param gsms: list -- list of GSM objects
        :param database: GEODatabase -- Database from SOFT file
        """

        gpls = {} if gpls is None else gpls
        if not isinstance(gpls, dict):
            raise ValueError("GPLs should be a dictionary not a %s" % str(type(gpls)))
        gsms = {} if gsms is None else gsms
        if not isinstance(gsms, dict):
            raise ValueError("GSMs should be a dictionary not a %s" % str(type(gsms)))

        for gsm_name, gsm in iteritems(gsms):
            assert isinstance(gsm, GSM), "All GSMs should be of type GSM"
        for gpl_name, gpl in iteritems(gpls):
            assert isinstance(gpl, GPL), "All GPLs should be of type GPL"
        if database is not None:
            if not isinstance(database, GEODatabase):
                raise ValueError("Database should be a GEODatabase not a %s" % str(type(database)))

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
                            char_type, char_value = [char[0], ": ".join(char[1:])]
                            tmp[key + "." + str(i) + "." + char_type] = char_value
                    else:
                        tmp[key] = ",".join(value)
                pheno_data[gsm_name] = tmp
            self._phenotype_data = DataFrame(pheno_data).T
        return self._phenotype_data

    def merge_and_average(self, platform, expression_column, group_by_column,
                          force=False, merge_on_column=None, gsm_on=None, gpl_on=None):
        """For given platform prepare the DataFrame with all the samples present in the GSE
        annotated with given column from platform and averaged over this column.

        :param platform: str or GPL -- GPL platform to use
        :param expression_column: str -- column name in which "expressions" are represented
        :param group_by_column: str -- the data will be grouped and averaged over this column and only this column will be kept
        :param force: bool -- if the name of the GPL does not match the platform name in GSM proceed anyway
        :param merge_on_column: str -- column to merge the data on - should be present in both GSM and GPL
        :param gsm_on: str -- in the case columns to merge are different in GSM and GPL use this column in GSM
        :param gpl_on: str -- in the case columns to merge are different in GSM and GPL use this column in GPL
        :returns: pandas.DataFrame

        """

        if isinstance(platform, str):
            gpl = self.gpls[platform]
        elif isinstance(platform, GPL):
            gpl = platform
        else:
            raise ValueError("Platform has to be of type GPL or string with key for platform in GSE")

        data = []
        for gsm in self.gsms.values():
            if gpl.name == gsm.metadata['platform_id'][0]:
                data.append(gsm.annotate_and_average(gpl=gpl,
                                                     merge_on_column=merge_on_column,
                                                     expression_column=expression_column,
                                                     group_by_column=group_by_column,
                                                     force=force,
                                                     gpl_on=gpl_on,
                                                     gsm_on=gsm_on))
        if len(data) == 0:
            stderr.write("No samples for the platform were found\n")
            return None
        elif len(data) == 1:
            return data[0]
        else:
            return data[0].join(data[1:])

    def pivot_samples(self, values, index="ID_REF"):
        """Construct a table in which columns (names) are the samples, index
        is a specified column eg. ID_REF and values in the columns are of one
        specified type.

        :param values: str -- column name present in the GSMs (all)
        :param index: str -- column name that will become an index in pivoted table
        :returns: pandas.DataFrame

        """
        data = []
        for gsm in self.gsms.values():
            tmp_data = gsm.table.copy()
            tmp_data["name"] = gsm.name
            data.append(tmp_data)
        ndf = concat(data).pivot(index=index, values=values, columns="name")
        return ndf

    def pivot_and_annotate(self, values, gpl, annotation_column, gpl_on="ID", gsm_on="ID_REF"):
        """Annotate GSM with provided GPL

        :param gpl: pandas.DataFrame or GPL -- a Platform or DataFrame to annotate with
        :param annotation_column: str -- column in table for annotation
        :param gsm_on: str -- use this column in GSM to merge, defaults to ID_REF
        :param gpl_on: str -- use this column in GPL to merge, defaults to ID
        :returns: pandas.DataFrame

        """
        if isinstance(gpl, GPL):
            annotation_table = gpl.table
        elif isinstance(gpl, DataFrame):
            annotation_table = gpl
        else:
            raise TypeError("gpl should be a GPL object or a pandas.DataFrame")
        pivoted_samples = self.pivot_samples(values=values, index=gsm_on)
        ndf = pivoted_samples.reset_index().merge(annotation_table[[gpl_on, annotation_column]],
                                                  left_on=gsm_on,
                                                  right_on=gpl_on).set_index(gsm_on)
        del ndf[gpl_on]
        ndf.columns.name = 'name'
        return ndf

    def download_supplementary_files(self, directory='series', download_sra=True, sra_filetype='fasta', email=None):
        """@todo: Docstring for download_supplementary_files.

        :param directory: directory to download the data (in this directory function will create
                          new directory with the files), by default this will be named with the series name + _Supp
        :param download_sra: bool - indicates whether to download SRA raw data too, defaults to True
        :param sra_filetype: indicates what file type to download if we specified SRA, can be sra, fasta or fastq
        :param email: e-mail that will be provided to the Entrez, defaults to None
        """
        if directory == 'series':
            dirpath = os.path.abspath(self.get_accession() + "_Supp")
            utils.mkdir_p(dirpath)
        else:
            dirpath = os.path.abspath(directory)
            utils.mkdir_p(dirpath)
        for gsmname, gsm in iteritems(self.gsms):
            gsm.download_supplementary_files(email=email, download_sra=download_sra, sra_filetype=sra_filetype, directory=dirpath)

    def download_SRA(self,  email, directory='series', filetype='sra', filterby=None):
        """Download SRA files for each GSM in series

        :param email: e-mail that will be provided to the Entrez
        :param directory: directory to save the data (defaults to the 'series' which saves the data to the
                          directory with the name of the series + '_SRA' ending)
        :param filetype: can be sra, fasta, or fastq - for fasta or fastq SRA-Toolkit need to be installed
        :param filterby: filter GSM objects, argument is a function that operates on GSM object  and return bool eg. lambda x: "brain" not in x.name

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

        for gsm in gsms_to_use:
            stderr.write("Downloading %s files for %s series\n" % (filetype, gsm.name))
            gsm.download_SRA(email=email, filetype=filetype, directory=dirpath)

    def _get_object_as_soft(self):
        """
         Return object as SOFT formated string.
        """
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
        return str("<%s: %s - %i SAMPLES, %i PLATFORM(s)>" % (self.geotype, self.name, len(self.gsms), len(self.gpls)))

    def __repr__(self):
        return str("<%s: %s - %i SAMPLES, %i PLATFORM(s)>" % (self.geotype, self.name, len(self.gsms), len(self.gpls)))

