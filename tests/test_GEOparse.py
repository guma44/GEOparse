#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_GEOparse
----------------------------------

Tests for `GEOparse` module.
"""

import unittest
from os.path import abspath, dirname, isdir, isfile, join
from sys import path

import GEOparse as GEO
from GEOparse.GEOTypes import (
    GDS,
    GPL,
    GSE,
    GSM,
    DataIncompatibilityException,
    GDSSubset,
    GEODatabase,
    NoMetadataException,
)
from pandas import DataFrame, read_table
from pandas.testing import assert_frame_equal
from six import iteritems

download_geo = dirname(abspath(__file__))


class TestGSM(unittest.TestCase):
    """Test GSM class"""

    def setUp(self):
        self.header1 = ["a", "b", "c"]
        self.header2 = ["a", "b", "d"]
        self.columns_desc = [["first column"], ["second column"], ["third column"]]
        self.data = [[1, 2, 3], [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns_no_desc = DataFrame(self.columns_desc, self.header1)
        self.columns1 = DataFrame(
            self.columns_desc, self.header1, columns=["description"]
        )
        self.columns2 = DataFrame(
            self.columns_desc, self.header2, columns=["description"]
        )
        self.metadata = {"name": ["SAMPLE"]}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GSM(name="name", table=["a"], metadata=self.metadata, columns=self.columns1)
        with self.assertRaises(TypeError):
            GSM(name="name", table=self.table1, metadata=[], columns=self.columns1)
        with self.assertRaises(ValueError):
            GSM(name="name", table=self.table1, metadata=self.metadata, columns=[])
        with self.assertRaises(DataIncompatibilityException):
            GSM(
                name="name",
                table=self.table1,
                metadata=self.metadata,
                columns=self.columns2,
            )
        with self.assertRaises(ValueError):
            GSM(
                name="name",
                table=self.table1,
                metadata=self.metadata,
                columns=self.columns_no_desc,
            )
        GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )

    def test_simple_data(self):
        gsm = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        self.assertEqual(gsm.table.loc[0, "a"], 1)
        self.assertEqual(gsm.table.loc[1, "b"], 5)

    def test_get_geo_and_data(self):
        gsm = GEO.get_GEO(geo="GSM11805", destdir=download_geo)
        self.assertTrue(isinstance(gsm, GSM))
        self.assertEqual(gsm.get_accession(), "GSM11805")
        self.assertEqual(len(gsm.table.index), 22283)
        self.assertEqual(len(gsm.columns), 3)
        self.assertEqual(len(gsm.metadata.keys()), 28)

    def test_get_metadata_attribute(self):
        metadata = {"name": "SAMPLE", "samples": ["sam1", "sam2"]}
        gsm = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        gsm2 = GSM(
            name="name", table=self.table1, metadata=metadata, columns=self.columns1
        )
        with self.assertRaises(TypeError):
            gsm2.get_metadata_attribute("name")
        self.assertEqual(gsm.get_metadata_attribute("name"), "SAMPLE")
        self.assertEqual(gsm2.get_metadata_attribute("samples"), ["sam1", "sam2"])
        with self.assertRaises(NoMetadataException):
            gsm.get_metadata_attribute("dupa")

    def test_get_type(self):
        gsm = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        self.assertEqual(gsm.get_type(), None)

    def test_metadata_as_string(self):
        metadata = {"name": ["SAMPLE1"], "data": ["Normal"]}
        gsm = GSM(
            name="name", table=self.table1, metadata=metadata, columns=self.columns1
        )
        self.assertIn("!Sample_data = Normal", gsm._get_metadata_as_string())
        self.assertIn("!Sample_name = SAMPLE1", gsm._get_metadata_as_string())

    def test_get_table_as_string(self):
        gsm = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        table = (
            "!sample_table_begin\n"
            "a\tb\tc\n"
            "1\t2\t3\n"
            "4\t5\t6\n"
            "!sample_table_end"
        )
        self.assertEqual(gsm._get_table_as_string(), table)

    def test_get_columns_as_string(self):
        gsm = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        columns = "#a = first column\n" "#b = second column\n" "#c = third column"
        self.assertEqual(gsm._get_columns_as_string(), columns)

    def test_to_soft(self):
        gsm = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        soft = (
            "^SAMPLE = name\n"
            "!Sample_name = SAMPLE\n"
            "#a = first column\n"
            "#b = second column\n"
            "#c = third column\n"
            "!sample_table_begin\n"
            "a\tb\tc\n"
            "1\t2\t3\n"
            "4\t5\t6\n"
            "!sample_table_end"
        )
        self.assertEqual(gsm._get_object_as_soft(), soft)

    def test_annotate(self):
        gse = GEO.get_GEO(
            filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE"
        )
        gsm = gse.gsms["Triple-Fusion Transfected Embryonic Stem Cells Replicate 1"]
        result = read_table(join(download_geo, "test_gsm_annotated.tab"))
        gpl = gse.gpls[next(iter(gse.gpls))]
        assert_frame_equal(result, gsm.annotate(gpl, annotation_column="GB_ACC"))
        assert_frame_equal(result, gsm.annotate(gpl.table, annotation_column="GB_ACC"))
        with self.assertRaises(TypeError):
            gsm.annotate("platform", annotation_column="GB_ACC")
        gsm.annotate(gpl.table, annotation_column="GB_ACC", in_place=True)
        assert_frame_equal(result, gsm.table)

    def test_head(self):
        gsm = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        try:
            gsm.head()
        except Exception:
            self.fail("GSM.head() raised error!")

    def test_empty_line(self):
        try:
            GEO.get_GEO(filepath=join(download_geo, "GSM32878.txt"), geotype="GSM")
        except IndexError:
            self.fail("Empty line in the file causes an error.")

    def test_no_table(self):
        try:
            GEO.get_GEO(filepath=join(download_geo, "GSM2795971.txt"), geotype="GSM")
        except Exception:
            self.fail("No data in the file error.")

    def test_download_SRA(self):
        geo_id = "GSE63525"  # Hi-C dataset from Rao et al.

        def filterby(x):
            return (
                "HIC174" in x.metadata["title"][0] or "HIC173" in x.metadata["title"][0]
            )

        destdir = "./TMP_SOFT"
        gse = GEO.get_GEO(geo=geo_id, destdir=destdir)
        downloaded_paths = gse.download_SRA(
            "ljosudorrit@gmail.com",  # some unused e-mail, use your own for run
            directory=destdir,
            filetype="sra",
            filterby=filterby,
            silent=True,
            keep_sra=True,
        )
        self.assertTrue(isdir(destdir))
        self.assertEqual(len(downloaded_paths), 2)
        for k in downloaded_paths.keys():
            self.assertTrue(k in gse.gsms.keys())
        for k in ["GSM1551718", "GSM1551719"]:
            self.assertTrue(k in downloaded_paths.keys())
        for k in downloaded_paths.keys():
            for f in downloaded_paths[k]["SRA"]:
                self.assertTrue(isfile(f))

    def test_download_SRA_parallel_by_gsm(self):
        geo_id = "GSE63525"  # Hi-C dataset from Rao et al.

        def filterby(x):
            return (
                "HIC173" in x.metadata["title"][0]
                or "HIC174" in x.metadata["title"][0]
                or "HIC175" in x.metadata["title"][0]
            )

        destdir = "./TMP_SOFT_parallel_by_gsm"

        gse = GEO.get_GEO(geo=geo_id, destdir=destdir)
        gsms_to_use = [gsm for gsm in gse.gsms.values() if filterby(gsm)]
        downloaded_paths = dict()
        for gsm in gsms_to_use:
            downloaded_paths[gsm.name] = gsm.download_SRA(
                "ljosudorrit@gmail.com",  # some unused e-mail
                directory=destdir,
                nproc=3,
                return_list=False,
                filetype="sra",
                silent=True,
                keep_sra=True,
            )
        self.assertTrue(isdir(destdir))
        self.assertEqual(len(downloaded_paths), 3)
        for k in downloaded_paths.keys():
            self.assertTrue(k in gse.gsms.keys())
        for k in ["GSM1551718", "GSM1551719", "GSM1551720"]:
            self.assertTrue(k in downloaded_paths.keys())
        for k in downloaded_paths.keys():
            for f in downloaded_paths[k]["SRA"]:
                self.assertTrue(isfile(f))

    def test_download_SRA_parallel_by_sra(self):
        geo_id = "GSE63525"  # Hi-C dataset from Rao et al.

        def filterby(x):
            return (
                "HIC173" in x.metadata["title"][0]
                or "HIC174" in x.metadata["title"][0]
                or "HIC175" in x.metadata["title"][0]
            )

        destdir = "./TMP_SOFT_parallel_by_sra"
        gse = GEO.get_GEO(geo=geo_id, destdir=destdir)
        downloaded_paths = gse.download_SRA(
            "ljosudorrit@gmail.com",  # some unused e-mail
            directory=destdir,
            filetype="sra",
            filterby=filterby,
            silent=True,
            keep_sra=True,
            nproc=3,
        )
        print(downloaded_paths)
        self.assertTrue(isdir(destdir))
        self.assertEqual(len(downloaded_paths), 3)
        for k in downloaded_paths.keys():
            self.assertTrue(k in gse.gsms.keys())
        for k in ["GSM1551718", "GSM1551719", "GSM1551720"]:
            self.assertTrue(k in downloaded_paths.keys())
        for k in downloaded_paths.keys():
            for f in downloaded_paths[k]["SRA"]:
                self.assertTrue(isfile(f))


class TestGPL(unittest.TestCase):
    """Test GPL class"""

    def setUp(self):
        self.header1 = ["a", "b", "c"]
        self.header2 = ["a", "b", "d"]
        self.columns_desc = [["first column"], ["second column"], ["third column"]]
        self.data = [[1, 2, 3], [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(
            self.columns_desc, self.header1, columns=["description"]
        )
        self.columns2 = DataFrame(
            self.columns_desc, self.header2, columns=["description"]
        )
        self.metadata = {"name": "PLATFORM"}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GPL(name="name", table=["a"], metadata=self.metadata, columns=self.columns1)
        with self.assertRaises(TypeError):
            GPL(name="name", table=self.table1, metadata=[], columns=self.columns1)
        with self.assertRaises(ValueError):
            GPL(name="name", table=self.table1, metadata=self.metadata, columns=[])
        with self.assertRaises(DataIncompatibilityException):
            GPL(
                name="name",
                table=self.table1,
                metadata=self.metadata,
                columns=self.columns2,
            )
        GPL(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )

    def test_simple_data(self):
        gpl = GPL(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        self.assertEqual(gpl.table.loc[0, "a"], 1)
        self.assertEqual(gpl.table.loc[1, "b"], 5)

    def test_get_geo_and_data(self):
        gpl = GEO.get_GEO(geo="GPL96", destdir=download_geo)
        self.assertTrue(isinstance(gpl, GPL))
        self.assertEqual(gpl.name, "GPL96")
        self.assertEqual(gpl.get_accession(), "GPL96")
        self.assertEqual(len(gpl.table.index), 22283)
        self.assertEqual(len(gpl.columns), 16)

    def test_get_geo_gpl_sequencing(self):
        gpl = GEO.get_GEO(geo="GPL20082", destdir=download_geo, include_data=True)
        self.assertTrue(isinstance(gpl, GPL))
        self.assertEqual(gpl.get_accession(), "GPL20082")

        samples = [
            "GSM1662787",
            "GSM1662788",
            "GSM1662789",
            "GSM1662790",
            "GSM1662791",
            "GSM1677167",
            "GSM1859499",
            "GSM1875285",
        ]

        for sample in samples:
            self.assertTrue(sample in gpl.gsms)

        self.assertEqual(6, len(gpl.gses["GSE68087"].gsms))
        self.assertEqual(2, len(gpl.gses["GSE67974"].gsms))

    def test_get_geo_gpl_partially(self):
        partial = ["GSM1662787", "GSM1662789", "GSM1662791", "GSM1859499"]

        gpl = GEO.get_GEO(
            geo="GPL20082", destdir=download_geo, include_data=True, partial=partial
        )
        self.assertTrue(isinstance(gpl, GPL))
        self.assertEqual(gpl.get_accession(), "GPL20082")

        for gsm in gpl.gsms:
            self.assertTrue(gsm in partial)

        self.assertEqual(4, len(gpl.gsms))

    def test_get_geo_and_data_with_annotations(self):
        gpl = GEO.get_GEO(geo="GPL96", destdir=download_geo, annotate_gpl=True)
        self.assertTrue(isinstance(gpl, GPL))
        self.assertEqual(gpl.name, "GPL96")
        self.assertEqual(gpl.get_metadata_attribute("platform"), "GPL96")
        self.assertEqual(len(gpl.table.index), 22283)
        self.assertEqual(len(gpl.columns), 21)

    def test_duplicate_column(self):
        columns = [
            "ID",
            "COL",
            "ROW",
            "NAME",
            "SPOT_ID",
            "CONTROL_TYPE",
            "REFSEQ",
            "GB_ACC",
            "GENE",
            "GENE_SYMBOL",
            "GENE_NAME",
            "UNIGENE_ID",
            "ENSEMBL_ID",
            "TIGR_ID",
            "ACCESSION_STRING",
            "CHROMOSOMAL_LOCATION",
            "CYTOBAND",
            "DESCRIPTION",
            "GO_ID",
            "SEQUENCE",
            "SPOT_ID.1",
            "ORDER",
        ]
        columns2 = [
            "ID",
            "COL",
            "ROW",
            "NAME",
            "SPOT_ID",
            "CONTROL_TYPE",
            "ENSEMBL_ID",
            "GB_ACC",
            "GENE",
            "GENE_SYMBOL",
            "ENSEMBL_ID.1",
            "UNIGENE_ID",
            "ENSEMBL_ID.2",
            "TIGR_ID",
            "ACCESSION_STRING",
            "CHROMOSOMAL_LOCATION",
            "CYTOBAND",
            "DESCRIPTION",
            "GO_ID",
            "SEQUENCE",
            "SPOT_ID.1",
            "ORDER",
        ]
        gpl = GEO.get_GEO(filepath=join(download_geo, "GPL4133.txt"))
        self.assertEqual(list(gpl.columns.index), columns)
        gpl2 = GEO.get_GEO(filepath=join(download_geo, "GPL4134.txt"))
        self.assertEqual(list(gpl2.columns.index), columns2)

    def test_name(self):
        gpl = GEO.get_GEO(
            filepath=join(download_geo, "GPL20814_family.soft"), geotype="GPL"
        )
        self.assertEqual(gpl.name, "GPL20814")


class TestGDS(unittest.TestCase):
    """Test GDS class"""

    def setUp(self):
        self.header1 = ["a", "b", "c"]
        self.header2 = ["a", "b", "d"]
        self.columns_desc = [["first column"], ["second column"], ["third column"]]
        self.data = [[1, 2, 3], [4, 5, 6]]
        self.subsets = {"s1": GDSSubset(name="subset", metadata={"type": "subset"})}
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(
            self.columns_desc, self.header1, columns=["description"]
        )
        self.columns2 = DataFrame(
            self.columns_desc, self.header2, columns=["description"]
        )
        self.metadata = {"name": "DATASET"}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GDS(
                name="name",
                table=["a"],
                metadata=self.metadata,
                columns=self.columns1,
                subsets=self.subsets,
            )
        with self.assertRaises(TypeError):
            GDS(
                name="name",
                table=self.table1,
                metadata=[],
                columns=self.columns1,
                subsets=self.subsets,
            )
        with self.assertRaises(ValueError):
            GDS(
                name="name",
                table=self.table1,
                metadata=self.metadata,
                columns=[],
                subsets=self.subsets,
            )
        GDS(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
            subsets=self.subsets,
        )

    def test_simple_data(self):
        gsm = GDS(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
            subsets=self.subsets,
        )
        self.assertEqual(gsm.table.loc[0, "a"], 1)
        self.assertEqual(gsm.table.loc[1, "b"], 5)

    def test_get_geo_and_data(self):
        gds = GEO.get_GEO(geo="GDS507", destdir=download_geo)
        self.assertTrue(isinstance(gds, GDS))
        self.assertEqual(gds.name, "GDS507")
        self.assertEqual(len(gds.table.index), 22645)
        self.assertEqual(len(gds.table.columns), 19)
        self.assertEqual(
            len(gds.metadata.keys()), 16
        )  # we omit DATABASE and SUBSET ! entries
        self.assertEqual(len(gds.database.metadata.keys()), 5)
        for subset_name, subset in iteritems(gds.subsets):
            self.assertEqual(len(subset.metadata.keys()), 4)
            self.assertTrue(isinstance(subset, GDSSubset))


class TestGSE(unittest.TestCase):
    """Test GSE class"""

    def setUp(self):
        self.header1 = ["a", "b", "c"]
        self.header2 = ["a", "b", "d"]
        self.columns_desc = [["first column"], ["second column"], ["third column"]]
        self.data = [[1, 2, 3], [4, 5, 6]]
        self.table1 = DataFrame(self.data, columns=self.header1)
        self.table2 = DataFrame(self.data, columns=self.header2)
        self.columns1 = DataFrame(
            self.columns_desc, self.header1, columns=["description"]
        )
        self.columns2 = DataFrame(
            self.columns_desc, self.header2, columns=["description"]
        )
        self.metadata = {"name": "PLATFORM"}
        self.gpl = GPL(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        self.gsm1 = GSM(
            name="name",
            table=self.table1,
            metadata=self.metadata,
            columns=self.columns1,
        )
        self.gsm2 = GSM(
            name="name",
            table=self.table2,
            metadata=self.metadata,
            columns=self.columns2,
        )
        self.gsms = {"a": self.gsm1, "b": self.gsm2}
        self.gpls = {"a": self.gpl}

    def test_creation_of_object(self):
        with self.assertRaises(ValueError):
            GSE(name="name1", metadata=self.metadata, gpls=[], gsms=self.gsms)
        with self.assertRaises(ValueError):
            GSE(name="name2", metadata=self.metadata, gpls=self.gpls, gsms=[])
        with self.assertRaises(TypeError):
            GSE(name="name3", metadata=[], gpls=self.gpls, gsms=self.gsms)
        GSE(name="name4", metadata=self.metadata, gpls=self.gpls, gsms=self.gsms)

    def test_soft_format_gse(self):
        print(download_geo)
        gse = GEO.get_GEO(geo="GSE1563", destdir=download_geo)
        self.assertTrue(isinstance(gse, GSE))
        self.assertEqual(gse.get_accession(), "GSE1563")
        self.assertEqual(len(gse.gsms.keys()), 62)
        self.assertEqual(len(gse.gpls.keys()), 1)
        self.assertEqual(len(gse.gpls[next(iter(gse.gpls))].table.index), 12625)
        self.assertEqual(len(gse.gsms[next(iter(gse.gsms))].table.index), 12625)
        for gsm_name, gsm in iteritems(gse.gsms):
            self.assertEqual(len(gsm.table.index), 12625)
            self.assertTrue(isinstance(gsm, GSM))
        for gpl_name, gpl in iteritems(gse.gpls):
            self.assertEqual(len(gpl.table.index), 12625)
            self.assertTrue(isinstance(gpl, GPL))

    def test_pivot_samples(self):
        gse = GEO.get_GEO(
            filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE"
        )
        result = read_table(
            join(download_geo, "test_sample_pivoted_by_value.tab"), index_col=0
        )
        result.columns.name = "name"
        assert_frame_equal(gse.pivot_samples("VALUE"), result)

    def test_merge_and_average(self):
        gse = GEO.get_GEO(
            filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE"
        )
        result = read_table(
            join(download_geo, "test_merged_by_id_and_averaged_by_gb_acc.tab"),
            index_col=0,
        )
        result = result.loc[
            sorted(result.index), sorted(result.columns)
        ]  # gse.gsms is a dict so the columns might be in different order
        merged = gse.merge_and_average(
            gse.gpls[next(iter(gse.gpls))],
            "VALUE",
            "GB_ACC",
            gpl_on="ID",
            gsm_on="ID_REF",
        )
        merged = merged[
            sorted(merged.columns)
        ]  # gse.gsms is a dict so the columns might be in different order
        assert_frame_equal(merged, result)
        with self.assertRaises(KeyError):
            gse.merge_and_average(
                "platform", "VALUE", "GB_ACC", gpl_on="ID", gsm_on="ID_REF"
            )
        with self.assertRaises(ValueError):
            gse.merge_and_average(
                ["platform"], "VALUE", "GB_ACC", gpl_on="ID", gsm_on="ID_REF"
            )

    def test_pivot_and_annotate(self):
        gse = GEO.get_GEO(
            filepath=join(download_geo, "soft_ex_family.txt"), geotype="GSE"
        )
        gpl = gse.gpls[next(iter(gse.gpls))]
        result = read_table(
            join(
                download_geo, "test_sample_pivoted_by_value_and_annotated_by_gbacc.tab"
            ),
            index_col=0,
        )
        result.columns.name = "name"
        pivoted = gse.pivot_and_annotate(
            values="VALUE", gpl=gpl, annotation_column="GB_ACC"
        )
        assert_frame_equal(result, pivoted)
        assert_frame_equal(
            gse.pivot_and_annotate(
                values="VALUE", gpl=gpl.table, annotation_column="GB_ACC"
            ),
            result,
        )
        with self.assertRaises(TypeError):
            gse.pivot_and_annotate(
                values="VALUE", gpl="gpl", annotation_column="GB_ACC"
            )

    def test_name(self):
        gse = GEO.get_GEO(
            filepath=join(download_geo, "GSE105845_family.soft"), geotype="GSE"
        )
        self.assertEqual(gse.name, "GSE105845")

    """
    TODO for GSE
    """

    def test_download_SRA(self):

        gse = GEO.get_GEO(geo="GSE1563", destdir=download_geo)
        self.assertTrue(isinstance(gse, GSE))
        self.assertEqual(gse.get_accession(), "GSE1563")
        self.assertEqual(len(gse.gsms.keys()), 62)
        self.assertEqual(len(gse.gpls.keys()), 1)
        self.assertEqual(len(gse.gpls[next(iter(gse.gpls))].table.index), 12625)
        self.assertEqual(len(gse.gsms[next(iter(gse.gsms))].table.index), 12625)
        for gsm_name, gsm in iteritems(gse.gsms):
            self.assertEqual(len(gsm.table.index), 12625)
            self.assertTrue(isinstance(gsm, GSM))
        for gpl_name, gpl in iteritems(gse.gpls):
            self.assertEqual(len(gpl.table.index), 12625)
            self.assertTrue(isinstance(gpl, GPL))


if __name__ == "__main__":
    unittest.main()
