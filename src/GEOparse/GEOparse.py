# -*- coding: utf-8 -*-
from collections import defaultdict
from itertools import groupby
from os import path
from re import match, split, sub
from tempfile import mkdtemp

from pandas import DataFrame, read_csv
from six import StringIO, iteritems

from . import utils
from .GEOTypes import GDS, GPL, GSE, GSM, GDSSubset, GEODatabase
from .logger import geoparse_logger as logger

try:
    from urllib.error import URLError
    from urllib.request import urlopen
except ImportError:
    from urllib2 import URLError, urlopen


class UnknownGEOTypeException(Exception):
    """Raised when the GEO type that do not correspond to any known."""

    pass


class NoEntriesException(Exception):
    """Raised when no entries could be found in the SOFT file."""

    pass


def get_GEO(
    geo=None,
    filepath=None,
    destdir="./",
    how="full",
    annotate_gpl=False,
    geotype=None,
    include_data=False,
    silent=False,
    aspera=False,
    partial=None,
    open_kwargs=None,
):
    """Get the GEO entry.

    The GEO entry is taken directly from the GEO database or read it from SOFT
    file.

    Args:
        geo (:obj:`str`): GEO database identifier.
        filepath (:obj:`str`): Path to local SOFT file. Defaults to None.
        destdir (:obj:`str`, optional): Directory to download data. Defaults to
            None.
        how (:obj:`str`, optional): GSM download mode. Defaults to "full".
            Possible options are: full, quick and brief
        annotate_gpl (:obj:`bool`, optional): Download the GPL annotation
            instead of regular GPL. If not available, fallback to regular GPL
            file. Defaults to False.
        geotype (:obj:`str`, optional): Type of GEO entry. By default it is
            inferred from the ID or the file name.
        include_data (:obj:`bool`, optional): Full download of GPLs including
            series and samples. Defaults to False.
        silent (:obj:`bool`, optional): Do not print anything. Defaults to
            False.
        aspera (:obj:`bool`, optional): EXPERIMENTAL Download using Aspera
            Connect. Follow Aspera instructions for further details. Defaults
            to False.
        partial (:obj:'iterable', optional): A list of accession IDs of GSMs
            to be partially extracted from GPL, works only if a file/accession
            is a GPL.
        open_kwargs (:obj:'dict', optional): A dict of kwargs that will be
            passed to `utils.smart_open` function.


    Returns:
        :obj:`GEOparse.BaseGEO`: A GEO object of given type.

    """
    if geo is None and filepath is None:
        raise Exception("You have to specify filename or GEO accession!")
    if geo is not None and filepath is not None:
        raise Exception("You can specify filename or GEO accession - not both!")
    if how not in {"full", "brief", "quick"}:
        raise Exception("Option 'how' can take only 'full', 'brief' or 'quick' values")

    if open_kwargs is None:
        open_kwargs = {}

    if silent:
        logger.setLevel(100)  # More than critical

    if filepath is None:
        filepath, geotype = get_GEO_file(
            geo,
            destdir=destdir,
            how=how,
            annotate_gpl=annotate_gpl,
            include_data=include_data,
            silent=silent,
            aspera=aspera,
        )
    else:
        if geotype is None:
            geotype = path.basename(filepath)[:3]

    logger.info("Parsing %s: " % filepath)
    if geotype.upper() == "GSM":
        return parse_GSM(filepath, open_kwargs=open_kwargs)
    elif geotype.upper() == "GSE":
        return parse_GSE(filepath, open_kwargs=open_kwargs)
    elif geotype.upper() == "GPL":
        return parse_GPL(filepath, partial=partial, open_kwargs=open_kwargs)
    elif geotype.upper() == "GDS":
        return parse_GDS(filepath, open_kwargs=open_kwargs)
    else:
        raise ValueError(
            ("Unknown GEO type: %s. Available types: GSM, GSE, " "GPL and GDS.")
            % geotype.upper()
        )


def get_GEO_file(
    geo,
    destdir=None,
    annotate_gpl=False,
    how="full",
    include_data=False,
    silent=False,
    aspera=False,
):
    """Download corresponding SOFT file given GEO accession.

    Args:
        geo (:obj:`str`): GEO database identifier.
        destdir (:obj:`str`, optional): Directory to download data. Defaults to
            None.
        annotate_gpl (:obj:`bool`, optional): Download the GPL annotation
            instead of regular GPL. If not available, fallback to regular GPL
            file. Defaults to False.
        how (:obj:`str`, optional): GSM download mode. Defaults to "full".
        include_data (:obj:`bool`, optional): Full download of GPLs including
            series and samples. Defaults to False.
        silent (:obj:`bool`, optional): Do not print anything. Defaults to
            False.
        aspera (:obj:`bool`, optional): EXPERIMENTAL Download using Aspera
            Connect. Follow Aspera instructions for further details. Defaults
            to False.

    Returns:
        :obj:`2-tuple` of :obj:`str` and :obj:`str`: Path to downloaded file and
        and the type of GEO object.

    """
    geo = geo.upper()
    geotype = geo[:3]
    range_subdir = sub(r"\d{1,3}$", "nnn", geo)
    if destdir is None:
        tmpdir = mkdtemp()
        logger.info(
            "No destination directory specified."
            " Temporary files will be downloaded at %s" % tmpdir
        )
    else:
        tmpdir = destdir
        utils.mkdir_p(tmpdir)

    if geotype == "GDS":
        gseurl = (
            "ftp://ftp.ncbi.nlm.nih.gov/geo/"
            "{root}/{range_subdir}/{record}/soft/{record_file}"
        )
        url = gseurl.format(
            root="datasets",
            range_subdir=range_subdir,
            record=geo,
            record_file="%s.soft.gz" % geo,
        )
        filepath = path.join(tmpdir, "{record}.soft.gz".format(record=geo))
    elif geotype == "GSE":
        if how == "full":
            gseurl = (
                "ftp://ftp.ncbi.nlm.nih.gov/geo/"
                "{root}/{range_subdir}/{record}/soft/{record_file}"
            )
            url = gseurl.format(
                root="series",
                range_subdir=range_subdir,
                record=geo,
                record_file="%s_family.soft.gz" % geo,
            )
            filepath = path.join(tmpdir, "{record}_family.soft.gz".format(record=geo))
        else:
            gseurl = (
                "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                "?targ=gsm&acc={record}&form=text&view={how}"
            )
            url = gseurl.format(record=geo, how=how)
            filepath = path.join(tmpdir, "{record}.txt".format(record=geo))
    elif geotype == "GSM":
        gsmurl = (
            "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
            "?targ=self&acc={record}&form=text&view={how}"
        )
        url = gsmurl.format(record=geo, how=how)
        filepath = path.join(tmpdir, "{record}.txt".format(record=geo))
    elif geotype == "GPL":
        if annotate_gpl:
            gplurl = (
                "ftp://ftp.ncbi.nlm.nih.gov/geo/"
                "{root}/{range_subdir}/{record}/annot/{record_file}"
            )
            url = gplurl.format(
                root="platforms",
                range_subdir=range_subdir,
                record=geo,
                record_file="%s.annot.gz" % geo,
            )
            filepath = path.join(tmpdir, "{record}.annot.gz".format(record=geo))
            if not path.isfile(filepath):
                try:
                    logger.info("Downloading %s to %s" % (url, filepath))
                    utils.download_from_url(url, filepath, silent=silent, aspera=aspera)
                    return filepath, geotype
                except URLError:
                    logger.info(
                        (
                            "Annotations for %s are not available, "
                            "trying submitter GPL"
                        )
                        % geo
                    )
            else:
                logger.info("File already exist: using local version.")
                return filepath, geotype

        if include_data:
            url = (
                "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/"
                "{0}/{1}/soft/{1}_family.soft.gz"
            ).format(range_subdir, geo)
            filepath = path.join(tmpdir, "{record}_family.soft.gz".format(record=geo))
        else:
            gplurl = (
                "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                "?targ=self&acc={record}&form=text&view={how}"
            )
            url = gplurl.format(record=geo, how=how)
            filepath = path.join(tmpdir, "{record}.txt".format(record=geo))
        if not path.isfile(filepath):
            logger.info("Downloading %s to %s" % (url, filepath))
            utils.download_from_url(url, filepath, silent=silent, aspera=aspera)
        else:
            logger.info("File already exist: using local version.")
        return filepath, geotype
    else:
        raise UnknownGEOTypeException("%s type is not known" % geotype)

    if not path.isfile(filepath):
        logger.info("Downloading %s to %s" % (url, filepath))
        utils.download_from_url(url, filepath, silent=silent, aspera=aspera)
    else:
        logger.info("File already exist: using local version.")

    return filepath, geotype


def __parse_entry(entry_line):
    """Parse the SOFT file entry name line that starts with '^', '!' or '#'.

    Args:
        entry_line (:obj:`str`): Line from SOFT  to be parsed.

    Returns:
        :obj:`2-tuple`: Type of entry, value of entry.

    """
    if entry_line.startswith("!"):
        entry_line = sub(r"!\w*?_", "", entry_line)
    else:
        entry_line = entry_line.strip()[1:]
    try:
        entry_type, entry_name = [i.strip() for i in entry_line.split("=", 1)]
    except ValueError:
        entry_type = [i.strip() for i in entry_line.split("=", 1)][0]
        entry_name = ""
    return entry_type, entry_name


def parse_entry_name(nameline):
    """Parse line that starts with ^ and assign the name to it.

    Args:
        nameline (:obj:`str`): A line to process.

    Returns:
        :obj:`str`: Entry name.

    """
    entry_type, entry_name = __parse_entry(nameline)
    return entry_name


def parse_metadata(lines):
    """Parse list of lines with metadata information from SOFT file.

    Args:
        lines (:obj:`Iterable`): Iterator over the lines.

    Returns:
        :obj:`dict`: Metadata from SOFT file.

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
    """Parse list of lines with columns description from SOFT file.

    Args:
        lines (:obj:`Iterable`): Iterator over the lines.

    Returns:
        :obj:`pandas.DataFrame`: Columns description.

    """
    data = []
    index = []
    for line in lines:
        line = line.rstrip()
        if line.startswith("#"):
            tmp = __parse_entry(line)
            data.append(tmp[1])
            index.append(tmp[0])

    return DataFrame(data, index=index, columns=["description"])


def parse_GDS_columns(lines, subsets):
    """Parse list of line with columns description from SOFT file of GDS.

    Args:
        lines (:obj:`Iterable`): Iterator over the lines.
        subsets (:obj:`dict` of :obj:`GEOparse.GDSSubset`): Subsets to use.

    Returns:
        :obj:`pandas.DataFrame`: Columns description.

    """
    data = []
    index = []
    for line in lines:
        line = line.rstrip()
        if line.startswith("#"):
            tmp = __parse_entry(line)
            data.append(tmp[1])
            index.append(tmp[0])

    df = DataFrame(data, index=index, columns=["description"])
    subset_ids = defaultdict(dict)
    for subsetname, subset in iteritems(subsets):
        for expid in subset.metadata["sample_id"][0].split(","):
            try:
                subset_type = subset.get_type()
                subset_ids[subset_type][expid] = subset.metadata["description"][0]
            except Exception:
                logger.error(
                    "Error processing subsets: %s for subset %s"
                    % (subset.get_type(), subsetname)
                )

    return df.join(DataFrame(subset_ids))


def parse_table_data(lines):
    """ "Parse list of lines from SOFT file into DataFrame.

    Args:
        lines (:obj:`Iterable`): Iterator over the lines.

    Returns:
        :obj:`pandas.DataFrame`: Table data.

    """
    # filter lines that do not start with symbols
    data = "\n".join(
        [i.rstrip() for i in lines if not i.startswith(("^", "!", "#")) and i.rstrip()]
    )
    if data:
        return read_csv(StringIO(data), index_col=None, sep="\t")
    else:
        return DataFrame()


def parse_GSM(filepath, entry_name=None, open_kwargs=None):
    """Parse GSM entry from SOFT file.

    Args:
        filepath (:obj:`str` or :obj:`Iterable`): Path to file with 1 GSM entry
            or list of lines representing GSM from GSE file.
        entry_name (:obj:`str`, optional): Name of the entry. By default it is
            inferred from the data.
        open_kwargs (:obj:'dict', optional): A dict of kwargs that will be
            passed to `utils.smart_open` function.

    Returns:
        :obj:`GEOparse.GSM`: A GSM object.

    """
    if open_kwargs is None:
        open_kwargs = {}
    if isinstance(filepath, str):
        with utils.smart_open(filepath, **open_kwargs) as f:
            soft = []
            has_table = False
            for line in f:
                if "_table_begin" in line or (not line.startswith(("^", "!", "#"))):
                    has_table = True
                soft.append(line.rstrip())
    else:
        soft = []
        has_table = False
        for line in filepath:
            if "_table_begin" in line or (not line.startswith(("^", "!", "#"))):
                has_table = True
            soft.append(line.rstrip())

    if entry_name is None:
        sets = [i for i in soft if i.startswith("^")]
        if len(sets) > 1:
            raise Exception("More than one entry in GPL")
        if len(sets) == 0:
            raise NoEntriesException(
                "No entries found. Check the if accession is correct!"
            )
        entry_name = parse_entry_name(sets[0])

    columns = parse_columns(soft)
    metadata = parse_metadata(soft)
    if has_table:
        table_data = parse_table_data(soft)
    else:
        table_data = DataFrame()

    gsm = GSM(name=entry_name, table=table_data, metadata=metadata, columns=columns)

    return gsm


def parse_GPL(filepath, entry_name=None, partial=None, open_kwargs=None):
    """Parse GPL entry from SOFT file.

    Args:
        filepath (:obj:`str` or :obj:`Iterable`): Path to file with 1 GPL entry
            or list of lines representing GPL from GSE file.
        entry_name (:obj:`str`, optional): Name of the entry. By default it is
            inferred from the data.
        partial (:obj:'iterable', optional): A list of accession IDs of GSMs
            to be partially extracted from GPL, works only if a file/accession
            is a GPL.
        open_kwargs (:obj:'dict', optional): A dict of kwargs that will be
            passed to `utils.smart_open` function.

    Returns:
        :obj:`GEOparse.GPL`: A GPL object.

    """
    gsms = {}
    gses = {}
    gpl_soft = []
    has_table = False
    gpl_name = entry_name
    database = None
    if open_kwargs is None:
        open_kwargs = {}
    if isinstance(filepath, str):
        with utils.smart_open(filepath, **open_kwargs) as soft:
            groupper = groupby(soft, lambda x: x.startswith("^"))
            for is_new_entry, group in groupper:
                if is_new_entry:
                    entry_type, entry_name = __parse_entry(next(group))
                    logger.debug("%s: %s" % (entry_type.upper(), entry_name))
                    if entry_type == "SERIES":
                        is_data, data_group = next(groupper)
                        gse_metadata = parse_metadata(data_group)

                        gses[entry_name] = GSE(name=entry_name, metadata=gse_metadata)
                    elif entry_type == "SAMPLE":
                        if partial and entry_name not in partial:
                            continue
                        is_data, data_group = next(groupper)
                        gsms[entry_name] = parse_GSM(data_group, entry_name)
                    elif entry_type == "DATABASE":
                        is_data, data_group = next(groupper)
                        database_metadata = parse_metadata(data_group)
                        database = GEODatabase(
                            name=entry_name, metadata=database_metadata
                        )

                    elif entry_type == "PLATFORM" or entry_type == "Annotation":
                        gpl_name = entry_name
                        is_data, data_group = next(groupper)
                        has_gpl_name = gpl_name or gpl_name is None
                        for line in data_group:
                            if "_table_begin" in line or not line.startswith(
                                ("^", "!", "#")
                            ):
                                has_table = True
                            if not has_gpl_name:
                                if match(r"!Annotation_platform\s*=\s*", line):
                                    gpl_name = split(r"\s*=\s*", line)[-1].strip()
                                    has_gpl_name = True
                            gpl_soft.append(line)
                    else:
                        raise RuntimeError(
                            "Cannot parse {etype}. Unknown for GPL.".format(
                                etype=entry_type
                            )
                        )
    else:
        for line in filepath:
            if "_table_begin" in line or (not line.startswith(("^", "!", "#"))):
                has_table = True
            gpl_soft.append(line.rstrip())

    columns = None
    try:
        columns = parse_columns(gpl_soft)
    except Exception:
        pass
    metadata = parse_metadata(gpl_soft)

    if has_table:
        table_data = parse_table_data(gpl_soft)
    else:
        table_data = DataFrame()

    gpl = GPL(
        name=gpl_name,
        gses=gses,
        gsms=gsms,
        table=table_data,
        metadata=metadata,
        columns=columns,
        database=database,
    )

    # link samples to series, if these were present in the GPL soft file
    for gse_id, gse in gpl.gses.items():
        for gsm_id in gse.metadata.get("sample_id", []):
            if gsm_id in gpl.gsms:
                gpl.gses[gse_id].gsms[gsm_id] = gpl.gsms[gsm_id]

    return gpl


def parse_GSE(filepath, open_kwargs=None):
    """Parse GSE SOFT file.

    Args:
        filepath (:obj:`str`): Path to GSE SOFT file.
        open_kwargs (:obj:'dict', optional): A dict of kwargs that will be
            passed to `utils.smart_open` function.

    Returns:
        :obj:`GEOparse.GSE`: A GSE object.

    """
    gpls = {}
    gsms = {}
    series_counter = 0
    database = None
    metadata = {}
    gse_name = None
    with utils.smart_open(filepath, **open_kwargs) as soft:
        groupper = groupby(soft, lambda x: x.startswith("^"))
        for is_new_entry, group in groupper:
            if is_new_entry:
                entry_type, entry_name = __parse_entry(next(group))
                logger.debug("%s: %s" % (entry_type.upper(), entry_name))
                if entry_type == "SERIES":
                    gse_name = entry_name
                    series_counter += 1
                    if series_counter > 1:
                        raise Exception(
                            "GSE file should contain only one series entry!"
                        )
                    is_data, data_group = next(groupper)
                    message = (
                        "The key is not False, probably there is an "
                        "error in the SOFT file"
                    )
                    assert not is_data, message
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
                    logger.error("Cannot recognize type %s" % entry_type)
    gse = GSE(name=gse_name, metadata=metadata, gpls=gpls, gsms=gsms, database=database)
    return gse


def parse_GDS(filepath, open_kwargs=None):
    """Parse GDS SOFT file.

    Args:
        filepath (:obj:`str`): Path to GDS SOFT file.
        open_kwargs (:obj:'dict', optional): A dict of kwargs that will be
            passed to `utils.smart_open` function.

    Returns:
        :obj:`GEOparse.GDS`: A GDS object.

    """
    dataset_lines = []
    subsets = {}
    database = None
    dataset_name = None
    with utils.smart_open(filepath, **open_kwargs) as soft:
        groupper = groupby(soft, lambda x: x.startswith("^"))
        for is_new_entry, group in groupper:
            if is_new_entry:
                entry_type, entry_name = __parse_entry(next(group))
                logger.debug("%s: %s" % (entry_type.upper(), entry_name))
                if entry_type == "SUBSET":
                    is_data, data_group = next(groupper)
                    message = (
                        "The key is not False, probably there is an "
                        "error in the SOFT file"
                    )
                    assert not is_data, message
                    subset_metadata = parse_metadata(data_group)
                    subsets[entry_name] = GDSSubset(
                        name=entry_name, metadata=subset_metadata
                    )
                elif entry_type == "DATABASE":

                    is_data, data_group = next(groupper)
                    message = (
                        "The key is not False, probably there is an "
                        "error in the SOFT file"
                    )
                    assert not is_data, message
                    database_metadata = parse_metadata(data_group)
                    database = GEODatabase(name=entry_name, metadata=database_metadata)
                elif entry_type == "DATASET":
                    is_data, data_group = next(groupper)
                    dataset_name = entry_name
                    for line in data_group:
                        dataset_lines.append(line.rstrip())
                else:
                    logger.error("Cannot recognize type %s" % entry_type)

    metadata = parse_metadata(dataset_lines)
    columns = parse_GDS_columns(dataset_lines, subsets)
    table = parse_table_data(dataset_lines)
    return GDS(
        name=dataset_name,
        metadata=metadata,
        columns=columns,
        table=table,
        subsets=subsets,
        database=database,
    )
