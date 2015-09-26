"""
Classes that represent different GEO entities
"""

from pandas import DataFrame, concat
from sys import stderr, stdout
import abc
import gzip


class DataIncompatibilityException(Exception): pass
class NoMetadataException(Exception): pass


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
        for metaname, meta in self.metadata.iteritems():
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
        print self._get_metadata_as_string()

    def to_soft(self, path_or_handle, as_gzip=False):
        """Save the object in a SOFT format.

        :param path_or_handle: path or handle to output file
        :param gzip: bool -- save as gzip
        """
        if isinstance(path_or_handle, str):
            if as_gzip:
                with gzip.open(path_or_handle, 'wb') as outfile:
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
        print self.columns

    def show_table(self, number_of_lines=5):
        """
        Show few lines of the table the table as pandas.DataFrame

        :param number_of_lines: int -- number of lines to show
        """
        print self.table.head(number_of_lines)

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


class GPL(SimpleGEO):

    """Class that represents platform from GEO database"""

    geotype = "PLATFORM"


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

        for subset_name, subset in subsets.iteritems():
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

    def __init__(self, name, metadata, gpls, gsms, database=None):
        """Initialize GSE

        :param name: str -- name of the object
        :param metadata: dict -- metadata information
        :param gpls: list -- list of GPL objects
        :param gsms: list -- list of GSM objects
        :param database: GEODatabase -- Database from SOFT file
        """

        if not isinstance(gpls, dict):
            raise ValueError("GPLs should be a dictionary not a %s" % str(type(gpls)))
        if not isinstance(gsms, dict):
            raise ValueError("GSMs should be a dictionary not a %s" % str(type(gsms)))

        for gsm_name, gsm in gsms.iteritems():
            assert isinstance(gsm, GSM), "All GSMs should be of type GSM"
        for gpl_name, gpl in gpls.iteritems():
            assert isinstance(gpl, GPL), "All GPLs should be of type GPL"
        if database is not None:
            if not isinstance(database, GEODatabase):
                raise ValueError("Database should be a GEODatabase not a %s" % str(type(database)))

        BaseGEO.__init__(self, name=name, metadata=metadata)

        self.gpls = gpls
        self.gsms = gsms
        self.database = database

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
        return ndf

    def _get_object_as_soft(self):
        """
         Return object as SOFT formated string.
        """
        soft = []
        if self.database is not None:
            soft.append(self.database._get_object_as_soft())
        soft += ["^%s = %s" % (self.geotype, self.name),
                 self._get_metadata_as_string()]
        for gsm in self.gsms.itervalues():
            soft.append(gsm._get_object_as_soft())
        for gpl in self.gpls.itervalues():
            soft.append(gpl._get_object_as_soft())

        return "\n".join(soft)

    def __str__(self):
        return str("<%s: %s - %i SERIES, %i PLATFORM(s)>" % (self.geotype, self.name, len(self.gsms), len(self.gpls)))

    def __repr__(self):
        return str("<%s: %s - %i SERIES, %i PLATFORM(s)>" % (self.geotype, self.name, len(self.gsms), len(self.gpls)))
