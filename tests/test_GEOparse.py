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
from pandas import DataFrame

path.append(join(abspath(__file__), ".."))
from GEOparse.GEOTypes import GSE, GSM, GPL, GDS, DataIncompatibilityException, GDSSubset, GEODatabase
import GEOparse as GEO

download_geo = dirname(abspath(__file__))

class TestGSM(unittest.TestCase):

    """Test GSM class"""

    def setUp(self):
        self.header1 = ['a', 'b', 'c']
        self.header2 = ['a', 'b', 'd']
        self.columns_desc = [['first column', 'second column', 'third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(self.columns_desc, self.header1)
        self.columns2 = DataFrame(self.columns_desc, self.header2)
        self.metadata = {'name': 'SAMPLE'}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GSM(name='name', table=['a'], metadata=self.metadata, columns=self.columns1)
        with self.assertRaises(ValueError):
            GSM(name='name', table=self.table1, metadata=[], columns=self.columns1)
        with self.assertRaises(ValueError):
            GSM(name='name', table=self.table1, metadata=self.metadata, columns=[])
        with self.assertRaises(DataIncompatibilityException):
            GSM(name='name', table=self.table1, metadata=self.metadata, columns=self.columns2)
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

class TestGPL(unittest.TestCase):

    """Test GPL class"""

    def setUp(self):
        self.header1 = ['a', 'b', 'c']
        self.header2 = ['a', 'b', 'd']
        self.columns_desc = [['first column', 'second column', 'third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(self.columns_desc, self.header1)
        self.columns2 = DataFrame(self.columns_desc, self.header2)
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

class TestGDS(unittest.TestCase):

    """Test GDS class"""

    def setUp(self):
        self.header1 = ['a', 'b', 'c']
        self.header2 = ['a', 'b', 'd']
        self.columns_desc = [['first column', 'second column', 'third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.subsets = {'s1': GDSSubset(name='subset', metadata={'type': 'subset'})}
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(self.columns_desc, self.header1)
        self.columns2 = DataFrame(self.columns_desc, self.header2)
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
        self.columns_desc = [['first column', 'second column', 'third column']]
        self.data = [[1, 2, 3],
                     [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(self.columns_desc, self.header1)
        self.columns2 = DataFrame(self.columns_desc, self.header2)
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

if __name__ == '__main__':
    unittest.main()

