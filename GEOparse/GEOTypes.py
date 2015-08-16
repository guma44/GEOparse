"""
Classes that represent different GEO entities
"""

from pandas import DataFrame
from sys import stderr

class DataIncompatibilityException(Exception): pass

class BaseGEO(object):

    def __init__(self, name, table, metadata, columns):
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
        if not isinstance(metadata, dict):
            raise ValueError("Metadata should be a dictionary not a %s" % str(type(metadata)))

        self.name = name
        self.table = table
        self.metadata = metadata
        self.columns = columns
        if self.columns.index.tolist() != self.table.columns.tolist():
            raise DataIncompatibilityException("Data columns do not match columns description index in %s" % (self.name))

    def get_accession(self):
        """Return accession ID of the sample
        :returns: str

        """
        return self.metadata["geo_accession"][0]


class GSM(BaseGEO):

    """Class that represents sample from GEO database"""

    def __init__(self, name, table, metadata, columns):
        """Initialize GSM sample

        :param name: str -- name of the object
        :param table: pandas.DataFrame -- table with the data from SOFT file
        :param metadata: dict -- metadata information
        :param columns: pandas.DataFrame -- description of the columns, number of columns, order, and names
        represented as index in this DataFrame has to be the same as table.columns.
        """

        BaseGEO.__init__(self, name=name, table=table, metadata=metadata, columns=columns)
        self.geotype = "SAMPLE"

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


class GPL(BaseGEO):

    """Class that represents platform from GEO database"""

    def __init__(self, name, table, metadata, columns):
        """Initialize GPL

        :param name: str -- name of the object
        :param table: pandas.DataFrame -- table with the data from SOFT file
        :param metadata: dict -- metadata information
        :param columns: pandas.DataFrame -- description of the columns, number of columns, order, and names
        represented as index in this DataFrame has to be the same as table.columns.
        """

        BaseGEO.__init__(self, name=name, table=table, metadata=metadata, columns=columns)
        self.geotype = "PLATFORM"


class GDS(BaseGEO):

    """Class that represents a dataset from GEO database"""

    def __init__(self, name, table, metadata, columns):
        """Initialize GDS

        :param name: str -- name of the object
        :param table: pandas.DataFrame -- table with the data from SOFT file
        :param metadata: dict -- metadata information
        :param columns: pandas.DataFrame -- description of the columns, number of columns, order, and names
        represented as index in this DataFrame has to be the same as table.columns.
        """

        BaseGEO.__init__(self, name=name, table=table, metadata=metadata, columns=columns)
        self.geotype = "DATASET"

class GSE(object):

    """Class representing GEO series"""

    def __init__(self, name, metadata, gpls, gsms):
        """Initialize GSE

        :param name: str -- name of the object
        :param metadata: dict -- metadata information
        :param gpls: list -- list of GPL objects
        :param gsms: list -- list of GSM objects

        """

        if not isinstance(metadata, dict):
            raise ValueError("Metadata should be a dictionary not a %s" % str(type(metadata)))
        if not isinstance(gpls, dict):
            raise ValueError("GPLs should be a dictionary not a %s" % str(type(gpls)))
        if not isinstance(gsms, dict):
            raise ValueError("GSMs should be a dictionary not a %s" % str(type(gsms)))

        for gsm_name, gsm in gsms.iteritems():
            assert isinstance(gsm, GSM), "All GSMs should be of type GSM"
        for gpl_name, gpl in gpls.iteritems():
            assert isinstance(gpl, GPL), "All GPLs should be of type GPL"

        self.name = name
        self.metadata = metadata
        self.gpls = gpls
        self.gsms = gsms

    def get_accession(self):
        """Return accession ID of the sample
        :returns: str

        """
        return self.metadata["geo_accession"][0]

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
