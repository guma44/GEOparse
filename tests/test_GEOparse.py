#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_GEOparse
----------------------------------

Tests for `GEOparse` module.
"""

from sys import path
from os.path import abspath, join, dirname
import unittest
from pandas import DataFrame, read_table
from pandas.util.testing import assert_frame_equal

path.append(join(abspath(__file__), ".."))
from GEOparse.GEOTypes import (GSE, GSM, GPL, GDS,
                               GDSSubset, GEODatabase,
                               DataIncompatibilityException,
                               NoMetadataException,
                               )
import GEOparse as GEO

download_geo = dirname(abspath(__file__))

class TestGSM(unittest.TestCase):

    """Test GSM class"""

    def setUp(self):
        self.header1 = ['a', 'b', 'c']
        self.header2 = ['a', 'b', 'd']
        self.columns_desc = [['first column'], ['second column'], ['third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns_no_desc = DataFrame(self.columns_desc, self.header1)
        self.columns1 = DataFrame(self.columns_desc, self.header1, columns=['description'])
        self.columns2 = DataFrame(self.columns_desc, self.header2, columns=['description'])
        self.metadata = {'name': ['SAMPLE']}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GSM(name='name', table=['a'], metadata=self.metadata, columns=self.columns1)
        with self.assertRaises(ValueError):
            GSM(name='name', table=self.table1, metadata=[], columns=self.columns1)
        with self.assertRaises(ValueError):
            GSM(name='name', table=self.table1, metadata=self.metadata, columns=[])
        with self.assertRaises(DataIncompatibilityException):
            GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns2)
        with self.assertRaises(ValueError):
            GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns_no_desc)
        GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)

    def test_simple_data(self):
        gsm = GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        self.assertEqual(gsm.table.ix[0, 'a'], 1)
        self.assertEqual(gsm.table.ix[1, 'b'], 5)

    def test_get_geo_and_data(self):
        gsm = GEO.get_GEO(geo="GSM11805", destdir=download_geo)
        self.assertTrue(isinstance(gsm, GSM))
        self.assertEqual(gsm.get_accession(), "GSM11805")
        self.assertEqual(len(gsm.table.index), 22283)
        self.assertEqual(len(gsm.columns), 3)
        self.assertEqual(len(gsm.metadata.keys()), 28)

    def test_get_metadata_attribute(self):
        metadata = {'name': 'SAMPLE', 'samples': ["sam1", "sam2"]}
        gsm = GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        gsm2 = GSM(name='name', table=self.table1, metadata=metadata, columns=self.columns1)
        with self.assertRaises(TypeError):
            gsm2.get_metadata_attribute('name')
        self.assertEqual(gsm.get_metadata_attribute('name'), "SAMPLE")
        self.assertEqual(gsm2.get_metadata_attribute('samples'), ["sam1", "sam2"])
        with self.assertRaises(NoMetadataException):
            gsm.get_metadata_attribute('dupa')

    def test_get_type(self):
        gsm = GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        self.assertEqual(gsm.get_type(), None)

    def test_metadata_as_string(self):
        metadata = {'name': ['SAMPLE1'], 'data': ["Normal"]}
        metadata_soft = ("!Sample_data = Normal\n"
                         "!Sample_name = SAMPLE1")
        gsm = GSM(name='name', table=self.table1, metadata=metadata, columns=self.columns1)
        self.assertEqual(gsm._get_metadata_as_string(), metadata_soft)

    def test_get_table_as_string(self):
        gsm = GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        table = ("!sample_table_begin\n"
                 "a\tb\tc\n"
                 "1\t2\t3\n"
                 "4\t5\t6\n"
                 "!sample_table_end")
        self.assertEqual(gsm._get_table_as_string(), table)

    def test_get_columns_as_string(self):
        gsm = GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        columns = ("#a = first column\n"
                   "#b = second column\n"
                   "#c = third column")
        self.assertEqual(gsm._get_columns_as_string(), columns)

    def test_to_soft(self):
        gsm = GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        soft = ("^SAMPLE = name\n"
                "!Sample_name = SAMPLE\n"
                "#a = first column\n"
                "#b = second column\n"
                "#c = third column\n"
                "!sample_table_begin\n"
                "a\tb\tc\n"
                "1\t2\t3\n"
                "4\t5\t6\n"
                "!sample_table_end")
        self.assertEqual(gsm._get_object_as_soft(), soft)

    def test_annotate(self):
        gse = GEO.get_GEO(filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE")
        gsm = gse.gsms["Triple-Fusion Transfected Embryonic Stem Cells Replicate 1"]
        result = read_table(join(download_geo, "test_gsm_annotated.tab"))
        gpl = gse.gpls[gse.gpls.keys()[0]]
        assert_frame_equal(result, gsm.annotate(gpl, annotation_column="GB_ACC"))
        assert_frame_equal(result, gsm.annotate(gpl.table, annotation_column="GB_ACC"))
        with self.assertRaises(TypeError):
            gsm.annotate("platform", annotation_column="GB_ACC")
        gsm.annotate(gpl.table, annotation_column="GB_ACC", in_place=True)
        assert_frame_equal(result, gsm.table)



class TestGPL(unittest.TestCase):

    """Test GPL class"""

    def setUp(self):
        self.header1 = ['a', 'b', 'c']
        self.header2 = ['a', 'b', 'd']
        self.columns_desc = [['first column'], ['second column'], ['third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(self.columns_desc, self.header1, columns=['description'])
        self.columns2 = DataFrame(self.columns_desc, self.header2, columns=['description'])
        self.metadata = {'name': 'PLATFORM'}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GPL(name='name', table=['a'], metadata=self.metadata, columns=self.columns1)
        with self.assertRaises(ValueError):
            GPL(name='name', table=self.table1, metadata=[], columns=self.columns1)
        with self.assertRaises(ValueError):
            GPL(name='name', table=self.table1, metadata=self.metadata, columns=[])
        with self.assertRaises(DataIncompatibilityException):
            GPL(name='name', table=self.table1, metadata=self.metadata, columns=self.columns2)
        GPL(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)

    def test_simple_data(self):
        gpl = GPL(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        self.assertEqual(gpl.table.ix[0, 'a'], 1)
        self.assertEqual(gpl.table.ix[1, 'b'], 5)

    def test_get_geo_and_data(self):
        gpl = GEO.get_GEO(geo="GPL96", destdir=download_geo)
        self.assertTrue(isinstance(gpl, GPL))
        self.assertEqual(gpl.get_accession(), "GPL96")
        self.assertEqual(len(gpl.table.index), 22283)
        self.assertEqual(len(gpl.columns), 16)

    def test_get_geo_and_data_with_annotations(self):
        gpl = GEO.get_GEO(geo="GPL96", destdir=download_geo, annotate_gpl=True)
        self.assertTrue(isinstance(gpl, GPL))
        self.assertEqual(gpl.get_metadata_attribute('platform'), "GPL96")
        self.assertEqual(len(gpl.table.index), 22283)
        self.assertEqual(len(gpl.columns), 21)

class TestGDS(unittest.TestCase):

    """Test GDS class"""

    def setUp(self):
        self.header1 = ['a', 'b', 'c']
        self.header2 = ['a', 'b', 'd']
        self.columns_desc = [['first column'], ['second column'], ['third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.subsets = {'s1': GDSSubset(name='subset', metadata={'type': 'subset'})}
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(self.columns_desc, self.header1, columns=['description'])
        self.columns2 = DataFrame(self.columns_desc, self.header2, columns=['description'])
        self.metadata = {'name': 'DATASET'}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GDS(name='name', table=['a'], metadata=self.metadata, columns=self.columns1, subsets=self.subsets)
        with self.assertRaises(ValueError):
            GDS(name='name', table=self.table1, metadata=[], columns=self.columns1, subsets=self.subsets)
        with self.assertRaises(ValueError):
            GDS(name='name', table=self.table1, metadata=self.metadata, columns=[], subsets=self.subsets)
        GDS(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1, subsets=self.subsets)

    def test_simple_data(self):
        gsm = GDS(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1, subsets=self.subsets)
        self.assertEqual(gsm.table.ix[0, 'a'], 1)
        self.assertEqual(gsm.table.ix[1, 'b'], 5)

    def test_get_geo_and_data(self):
        gds = GEO.get_GEO(geo="GDS507", destdir=download_geo)
        self.assertTrue(isinstance(gds, GDS))
        self.assertEqual(len(gds.table.index), 22645)
        self.assertEqual(len(gds.table.columns), 19)
        self.assertEqual(len(gds.metadata.keys()), 16) # we omit DATABASE and SUBSET ! entries
        self.assertEqual(len(gds.database.metadata.keys()), 5)
        for subset_name, subset in gds.subsets.iteritems():
            self.assertEqual(len(subset.metadata.keys()), 4)
            self.assertTrue(isinstance(subset, GDSSubset))

class TestGSE(unittest.TestCase):

    """Test GSE class"""

    def setUp(self):
        self.header1 = ['a', 'b', 'c']
        self.header2 = ['a', 'b', 'd']
        self.columns_desc = [['first column'], ['second column'], ['third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(self.columns_desc, self.header1, columns=['description'])
        self.columns2 = DataFrame(self.columns_desc, self.header2, columns=['description'])
        self.metadata = {'name': 'PLATFORM'}
        self.gpl = GPL(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        self.gsm1 = GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns1)
        self.gsm2 = GSM(name='name', table=self.table2, metadata=self.metadata, columns=self.columns2)
        self.gsms = {'a': self.gsm1,'b': self.gsm2}
        self.gpls = {'a': self.gpl}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GSE(name='name', metadata=self.metadata, gpls=[], gsms=self.gsms)
        with self.assertRaises(ValueError):
            GSE(name='name', metadata=self.metadata, gpls=self.gpls, gsms=[])
        with self.assertRaises(ValueError):
            GSE(name='name', metadata=[], gpls=self.gpls, gsms=self.gsms)
        GSE(name='name', metadata=self.metadata, gpls=self.gpls, gsms=self.gsms)

    def test_soft_format_gse(self):
        print download_geo
        gse = GEO.get_GEO(geo="GSE1563", destdir=download_geo)
        self.assertTrue(isinstance(gse, GSE))
        self.assertEqual(gse.get_accession(), "GSE1563")
        self.assertEqual(len(gse.gsms.keys()), 62)
        self.assertEqual(len(gse.gpls.keys()), 1)
        self.assertEqual(len(gse.gpls[gse.gpls.keys()[0]].table.index), 12625)
        self.assertEqual(len(gse.gsms[gse.gsms.keys()[0]].table.index), 12625)
        for gsm_name, gsm in gse.gsms.iteritems():
            self.assertEqual(len(gsm.table.index), 12625)
            self.assertTrue(isinstance(gsm, GSM))
        for gpl_name, gpl in gse.gpls.iteritems():
            self.assertEqual(len(gpl.table.index), 12625)
            self.assertTrue(isinstance(gpl, GPL))

    def test_pivot_samples(self):
        gse = GEO.get_GEO(filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE")
        result = read_table(join(download_geo, "test_sample_pivoted_by_value.tab"), index_col=0)
        result.columns.name = 'name'
        assert_frame_equal(gse.pivot_samples("VALUE"), result)

    def test_merge_and_average(self):
        gse = GEO.get_GEO(filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE")
        result = read_table(join(download_geo, "test_merged_by_id_and_averaged_by_gb_acc.tab"), index_col=0)
        result = result.ix[sorted(result.index), sorted(result.columns)]  # gse.gsms is a dict so the columns might be in different order
        merged = gse.merge_and_average(gse.gpls[gse.gpls.keys()[0]], "VALUE", "GB_ACC", gpl_on="ID", gsm_on="ID_REF")
        merged = merged[sorted(merged.columns)]  # gse.gsms is a dict so the columns might be in different order
        assert_frame_equal(merged, result)
        with self.assertRaises(KeyError):
            gse.merge_and_average("platform", "VALUE", "GB_ACC", gpl_on="ID", gsm_on="ID_REF")
        with self.assertRaises(ValueError):
            gse.merge_and_average(["platform"], "VALUE", "GB_ACC", gpl_on="ID", gsm_on="ID_REF")

    def test_pivot_and_annotate(self):
        gse = GEO.get_GEO(filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE")
        gpl = gse.gpls[gse.gpls.keys()[0]]
        result = read_table(join(download_geo, "test_sample_pivoted_by_value_and_annotated_by_gbacc.tab"), index_col=0)
        result.columns.name = 'name'
        pivoted = gse.pivot_and_annotate(values="VALUE", gpl=gpl, annotation_column="GB_ACC")
        assert_frame_equal(result, pivoted)
        assert_frame_equal(gse.pivot_and_annotate(values="VALUE", gpl=gpl.table, annotation_column="GB_ACC"),
                           result)
        with self.assertRaises(TypeError):
            gse.pivot_and_annotate(values="VALUE", gpl="gpl", annotation_column="GB_ACC")


if __name__ == '__main__':
    unittest.main()

