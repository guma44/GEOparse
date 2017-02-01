# -*- coding: utf-8 -*-
from os import path
from re import sub
from sys import stderr
from six import StringIO
from tempfile import mkdtemp
from itertools import groupby
from collections import defaultdict
from shutil import copyfileobj
try:
    from urllib.request import urlopen
    from urllib.error import URLError
except ImportError:
    from urllib2 import urlopen, URLError
from contextlib import closing
from pandas import DataFrame
import gzip
from .GEOTypes import GSE, GSM, GPL, GDS, GDSSubset, GEODatabase
from six import iteritems
import wgetter


class UnknownGEOTypeException(Exception):
    """
    Exception representing the GEO type that do not correspond to any known.
    """
    pass

class NoEntriesException(Exception):
    """
    Exception raised when no entries could be found in the SOFT file.
    """
    pass


def get_GEO(geo=None, filepath=None, destdir="./", how='full', annotate_gpl=False, geotype=None, include_data=False, silent=False):
    """Get the GEO entry directly from the GEO database or read it from SOFT file.

    :param geo: str -- GEO database identifier
    :param filepath: str -- path to local SOFT file
    :param destdir: str -- directory to download data
    :param how: str -- GSM download mode: full ...
    :param include_data: bool -- full download of GPLs including series and samples
    :param silent: bool -- don't print info
    :returns: GEOType object -- object according to specified GEO type

    """
    if (geo is None and filepath is None):
        raise Exception("You have to specify filename or GEO accession!")
    if (geo is not None and filepath is not None):
        raise Exception("You can specify filename or GEO accession - not both!")

    if filepath is None:
        filepath, geotype = get_GEO_file(geo, destdir=destdir, how=how, annotate_gpl=annotate_gpl, include_data=include_data)
    else:
        if geotype is None:
            geotype = filepath.split("/")[-1][:3]

    stderr.write("Parsing %s:\n" % filepath)
    if geotype.upper() == "GSM":
        return parse_GSM(filepath)
    elif geotype.upper() == "GSE":
        return parse_GSE(filepath)
    elif geotype.upper() == 'GPL':
        return parse_GPL(filepath, silent=silent)
    elif geotype.upper() == 'GDS':
        return parse_GDS(filepath)
    else:
        raise ValueError("Unknown GEO type: %s. Available types: GSM, GSE, GPL and GDS." % geotype.upper())


def get_GEO_file(geo, destdir=None, annotate_gpl=False, how="full",
        include_data=False):
    """Given GEO accession download corresponding SOFT file

    :param geo: str -- GEO database identifier
    :param destdir: str -- directory to download data
    :param how: str -- GSM download mode: full ...
    :param include_data: bool -- full download of GPLs including series and samples
    :returns: tuple -- path to downladed file, type of GEO object

    """
    geo = geo.upper()
    geotype = geo[:3]
    range_subdir = sub(r"\d{1,3}$", "nnn", geo)
    if destdir is None:
        tmpdir = mkdtemp()
        stderr.write("No destination directory specified."
                     " Temporary files will be downloaded at %s\n" % tmpdir)
    else:
        tmpdir = destdir

    if geotype == "GDS":
        gseurl = "ftp://ftp.ncbi.nlm.nih.gov/geo/{root}/{range_subdir}/{record}/soft/{record_file}"
        url = gseurl.format(root="datasets",
                            range_subdir=range_subdir,
                            record=geo,
                            record_file="%s.soft.gz" % geo)
        filepath = path.join(tmpdir, "{record}.soft.gz".format(record=geo))
    elif geotype == "GSE":
        gseurl = "ftp://ftp.ncbi.nlm.nih.gov/geo/{root}/{range_subdir}/{record}/soft/{record_file}"
        url = gseurl.format(root="series",
                            range_subdir=range_subdir,
                            record=geo,
                            record_file="%s_family.soft.gz" % geo)
        filepath = path.join(tmpdir, "{record}_family.soft.gz".format(record=geo))
    elif geotype == "GSM":
        gsmurl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={record}&form=text&view={how}"
        url = gsmurl.format(record=geo, how=how)
        filepath = path.join(tmpdir, "{record}.txt".format(record=geo))
    elif geotype == "GPL":
        if annotate_gpl:
            gplurl = "ftp://ftp.ncbi.nlm.nih.gov/geo/{root}/{range_subdir}/{record}/annot/{record_file}"
            url = gplurl.format(root="platforms",
                                range_subdir=range_subdir,
                                record=geo,
                                record_file="%s.annot.gz" % geo)
            filepath = path.join(tmpdir, "{record}.annot.gz".format(record=geo))
            if not path.isfile(filepath):
                try:
                    stderr.write("Downloading %s to %s\n" % (url, filepath))
                    fn = wgetter.download(url, outdir=path.dirname(filepath))
                    return filepath, geotype
                except URLError:
                    stderr.write("Annotations for %s are not available, trying submitter GPL\n" % geo)
            else:
                stderr.write("File already exist: using local version.\n")
                return filepath, geotype

        if include_data:
            url = "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/{0}/{1}/soft/{1}_family.soft.gz".format(
                    range_subdir,
                    geo)
            filepath = path.join(tmpdir, "{record}_family.soft.gz".format(record=geo))
        else:
            gplurl = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc={record}&form=text&view={how}"
            url = gplurl.format(record=geo, how=how)
            filepath = path.join(tmpdir, "{record}.soft".format(record=geo))
        if not path.isfile(filepath):
            stderr.write("Downloading %s to %s\n" % (url, filepath))
            fn = wgetter.download(url, outdir=path.dirname(filepath))
            stderr.write("\n")
        else:
            stderr.write("File already exist: using local version.\n")
        return filepath, geotype
    else:
        raise UnknownGEOTypeException("%s type is not known" % geotype)

    if not path.isfile(filepath):
        stderr.write("Downloading %s to %s\n" % (url, filepath))
        fn = wgetter.download(url, outdir=path.dirname(filepath))
        stderr.write("\n")
    else:
        stderr.write("File already exist: using local version.\n")

    return filepath, geotype


def __parse_entry(entry_line):
    """Parse the SOFT file entry name line that starts with '^', '!' or '#'.

    :param entry_line: str -- line from SOFT file
    :returns: tuple -- type, value

    """
    if entry_line.startswith("!"):
        entry_line = sub(r"!\w*?_", '', entry_line)
    else:
        entry_line = entry_line.strip()[1:]
    try:
        entry_type, entry_name = [i.strip() for i in entry_line.split("=", 1)]
    except ValueError:
        entry_type = [i.strip() for i in entry_line.split("=", 1)][0]
        entry_name = ''
    return entry_type, entry_name


def parse_entry_name(nameline):
    """Parse line that starts with ^ and assign the name to it

    :param nameline: str -- line to process
    :returns: str -- entry name

    """
    entry_type, entry_name = __parse_entry(nameline)
    return entry_name


def parse_metadata(lines):
    """Parse list of lines with metadata information from SOFT file

    :param lines: iterable -- iterator over lines
    :returns: dict -- metadata from SOFT file

    """
    meta = defaultdict(list)
    for line in lines:
        line = line.rstrip()
        if line.startswith("!"):
            if "_table_begin" in line or "_table_end" in line:
                continue
            key, value = __parse_entry(line)
            meta[key].append(value)

    return dict(meta)


def parse_columns(lines):
    """Parse list of line with columns description from SOFT file

    :param lines: iterable -- iterator over lines
    :returns: pandas.DataFrame -- columns description

    """
    data = []
    index = []
    for line in lines:
        line = line.rstrip()
        if line.startswith("#"):
            tmp = __parse_entry(line)
            data.append(tmp[1])
            index.append(tmp[0])

    return DataFrame(data, index=index, columns=['description'])


def parse_GDS_columns(lines, subsets):
    """Parse list of line with columns description from SOFT file
    of GDS (GEO Dataset)

    :param lines: iterable -- iterator over lines
    :returns: pandas.DataFrame -- columns description

    """
    data = []
    index = []
    for line in lines:
        line = line.rstrip()
        if line.startswith("#"):
            tmp = __parse_entry(line)
            data.append(tmp[1])
            index.append(tmp[0])

    df = DataFrame(data, index=index, columns=['description'])
    subset_ids = {"disease_state": {}, "individual": {}}
    for subsetname, subset in iteritems(subsets):
        for expid in subset.metadata["sample_id"][0].split(","):
            if subset.get_type() == "disease state":
                subset_ids["disease_state"][expid] = subset.metadata["description"][0]
            elif subset.get_type() == "individual":
                subset_ids["individual"][expid] = subset.metadata["description"][0]
            else:
                stderr.write("Unknown subset type: %s for subset %s\n" % (subset.get_type(), subsetname))

    return df.join(DataFrame(subset_ids))


def parse_table_data(lines):
    """Parse list of lines from SOFT file into DataFrame

    :param lines: iterable -- iterator over lines
    :returns: pandas.DataFrame -- table data

    """
    # filter lines that do not start with symbols
    data = "\n".join([i.rstrip() for i in lines if i[0] not in ("^", "!", "#")])
    return DataFrame.from_csv(StringIO(data), index_col=None, sep="\t")


def parse_GSM(filepath, entry_name=None):
    """Parse GSM entry from SOFT file

    :param filepath: str or iterable -- path to file with 1 GSM entry or list of lines representing
                                    GSM from GSE file
    :return: GSM object

    """
    if isinstance(filepath, str):
        if filepath[-2:] == "gz":
            mode = "rt"
            fopen = gzip.open
        else:
            mode = "r"
            fopen = open
        with fopen(filepath, mode) as f:
            soft = []
            has_table = False
            for line in f:
                if "_table_begin" in line or (line[0] not in ("^", "!", "#")):
                    has_table = True
                soft.append(line.rstrip())
    else:
        soft = []
        has_table = False
        for line in filepath:
            if "_table_begin" in line or (line[0] not in ("^", "!", "#")):
                has_table = True
            soft.append(line.rstrip())

    if entry_name is None:
        sets = [i for i in soft if i.startswith("^")]
        if len(sets) > 1:
            raise Exception("More than one entry in GPL")
        if len(sets) == 0:
            raise NoEntriesException("No entries found. Check the if accession is correct!")
        entry_name = parse_entry_name(sets[0])

    columns = parse_columns(soft)
    metadata = parse_metadata(soft)
    if has_table:
        table_data = parse_table_data(soft)
    else:
        table_data = DataFrame()

    gsm = GSM(name=entry_name,
              table=table_data,
              metadata=metadata,
              columns=columns)

    return gsm


def parse_GPL(filepath, entry_name=None, silent=False):
    """Parse GPL entry from SOFT file

    :param filepath: str or iterable -- path to file with 1 GPL entry or list of lines representing
                                    GPL from GSE file
    :param silent: bool -- don't print info
    :return: GPL object

    """
    gpls = {}
    gsms = {}
    gses = {}
    metadata = {}
    gpl_soft = []
    has_table = False
    if isinstance(filepath, str):
        if filepath[-2:] == "gz":
            mode = "rt"
            fopen = gzip.open
        else:
            mode = "r"
            fopen = open
        with fopen(filepath, mode) as soft:
            groupper = groupby(soft, lambda x: x.startswith("^"))
            for is_new_entry, group in groupper:
                if is_new_entry:
                    entry_type, entry_name = __parse_entry(next(group))
                    if not silent:
                        stderr.write(" - %s : %s\n" % (entry_type.upper(), entry_name))
                    if entry_type == "SERIES":
                        is_data, data_group = next(groupper)
                        gse_metadata = parse_metadata(data_group)

                        gses[entry_name] = GSE(name=entry_name, metadata=gse_metadata)
                    elif entry_type == "SAMPLE":
                        is_data, data_group = next(groupper)
                        gsms[entry_name] = parse_GSM(data_group, entry_name)
                    elif entry_type == "DATABASE":
                        is_data, data_group = next(groupper)
                        database_metadata = parse_metadata(data_group)
                        database = GEODatabase(name=entry_name, metadata=database_metadata)
                    else:
                        is_data, data_group = next(groupper)
                        for line in data_group:
                            if "_table_begin" in line or (line[0] not in ("^", "!", "#")):
                                has_table = True
                            gpl_soft.append(line)
    else:
         for line in filepath:
             if "_table_begin" in line or (line[0] not in ("^", "!", "#")):
                 has_table = True
             gpl_soft.append(line.rstrip())

    columns = None
    try:
        columns = parse_columns(gpl_soft)
    except:
        pass
    metadata = parse_metadata(gpl_soft)

    if has_table:
        table_data = parse_table_data(gpl_soft)
    else:
        table_data = DataFrame()

    gpl = GPL(name=entry_name,
              gses=gses,
              gsms=gsms,
              table=table_data,
              metadata=metadata,
              columns=columns,
              )

    # link samples to series, if these were present in the GPL soft file
    for gse_id,gse in gpl.gses.items():
        for gsm_id in gse.metadata.get("sample_id", []):
            if gsm_id in gpl.gsms:
                gpl.gses[gse_id].gsms[gsm_id] = gpl.gsms[gsm_id]

    return gpl

def parse_GSE(filepath):
    """Parse GSE from SOFT file

    :param filepath: str -- path to GSE SOFT file
    :return: GSE object
    """
    if filepath[-2:] == "gz":
        mode = "rt"
        fopen = gzip.open
    else:
        mode = "r"
        fopen = open
    gpls = {}
    gsms = {}
    series_counter = 0
    database = None
    metadata = {}
    with fopen(filepath, mode) as soft:
        groupper = groupby(soft, lambda x: x.startswith("^"))
        for is_new_entry, group in groupper:
            if is_new_entry:
                entry_type, entry_name = __parse_entry(next(group))
                stderr.write(" - %s : %s\n" % (entry_type.upper(), entry_name))
                if entry_type == "SERIES":
                    series_counter += 1
                    if series_counter > 1:
                        raise Exception("GSE file should contain only one series entry!")
                    is_data, data_group = next(groupper)
                    assert not is_data, "The key is not False, probably there is an error in the SOFT file"
                    metadata = parse_metadata(data_group)
                elif entry_type == "SAMPLE":
                    is_data, data_group = next(groupper)
                    gsms[entry_name] = parse_GSM(data_group, entry_name)
                elif entry_type == "PLATFORM":
                    is_data, data_group = next(groupper)
                    gpls[entry_name] = parse_GPL(data_group, entry_name)
                elif entry_type == "DATABASE":
                    is_data, data_group = next(groupper)
                    database_metadata = parse_metadata(data_group)
                    database = GEODatabase(name=entry_name, metadata=database_metadata)
                else:
                    stderr.write("Cannot recognize type %s\n" % entry_type)
    gse = GSE(name=entry_name,
              metadata=metadata,
              gpls=gpls,
              gsms=gsms,
              database=database)
    return gse


def parse_GDS(filepath):
    """
    Parse GDS from SOFT file

    :param gds: @todo
    :returns: @todo

    """
    if filepath[-2:] == "gz":
        mode = "rt"
        fopen = gzip.open
    else:
        mode = "r"
        fopen = open
    dataset_lines = []
    subsets = {}
    database = None
    with fopen(filepath, mode) as soft:
        groupper = groupby(soft, lambda x: x.startswith("^"))
        for is_new_entry, group in groupper:
            if is_new_entry:
                entry_type, entry_name = __parse_entry(next(group))
                stderr.write(" - %s : %s\n" % (entry_type.upper(), entry_name))
                if entry_type == "SUBSET":
                    is_data, data_group = next(groupper)
                    assert not is_data, "The key is not False, probably there is an error in the SOFT file"
                    subset_metadata = parse_metadata(data_group)
                    subsets[entry_name] = GDSSubset(name=entry_name, metadata=subset_metadata)
                elif entry_type == "DATABASE":

                    is_data, data_group = next(groupper)
                    assert not is_data, "The key is not False, probably there is an error in the SOFT file"
                    database_metadata = parse_metadata(data_group)
                    database = GEODatabase(name=entry_name, metadata=database_metadata)
                elif entry_type == "DATASET":
                    is_data, data_group = next(groupper)
                    dataset_name = entry_name
                    for line in data_group:
                        dataset_lines.append(line.rstrip())
                else:
                    stderr.write("Cannot recognize type %s\n" % entry_type)

    metadata = parse_metadata(dataset_lines)
    columns = parse_GDS_columns(dataset_lines, subsets)
    table = parse_table_data(dataset_lines)
    return GDS(name=dataset_name, metadata=metadata, columns=columns, table=table, subsets=subsets, database=database)
