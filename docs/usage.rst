========
Usage
========

To simplest example of usage::

    import GEOparse

    gse = GEOparse.get_GEO(geo="GSE1563", destdir="./")

    print()
    print("GSM example:")
    for gsm_name, gsm in gse.gsms.items():
        print("Name: ", gsm_name)
        print("Metadata:",)
        for key, value in gsm.metadata.items():
            print(" - %s : %s" % (key, ", ".join(value)))
        print ("Table data:",)
        print (gsm.table.head())
        break

    print()
    print("GPL example:")
    for gpl_name, gpl in gse.gpls.items():
        print("Name: ", gpl_name)
        print("Metadata:",)
        for key, value in gpl.metadata.items():
            print(" - %s : %s" % (key, ", ".join(value)))
        print("Table data:",)
        print(gpl.table.head())
        break

Working with GEO accession
^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to download and load to the memory the experiment with specific
GEO accession use :meth:`GEOparse.get_GEO`.

Eg. if you would like to download the *GSE1563* the code::

    import GEOparse
    gse = GEOparse.get_GEO(geo="GSE1563", destdir="./")

will check the GEO database for this accession ID and it will download it to specified
directory (by default this is CWD). The result will be loaded into :class:`GEOparse.GSE`.
For more information on how to work with GEOparse objects see :ref:`working-with-objects`.


Working with SOFT files
^^^^^^^^^^^^^^^^^^^^^^^

Reading SOFT file
-----------------

If you already downloaded some GDS, GSE, GSM or GPL you can read it into object by explicitly
specifying path to the file::

    import GEOparse
    gse = GEOparse.get_GEO(filepath="./GSE1563.soft.gz")

This will read the file into GSE object.

Crating and saving SOFT file
----------------------------

When doing own experiments one might want to submit the data to the GEO. In this case, one shuld
first create all necessary elements of given GEO object, create this object and use :func:`to_soft`
method to save it as SOFT files.

.. warning::

    As for now there is no implementation of any check procedure if the tables and metadata contain
    all necessary information required by GEO. This could be introduced in the future, however, now it is
    entirely up to the user to follow the rules of GEO.

Here is the example of a faked GSE that will be saved as SOFT. There is three main components of
GSE object: a dictionary of GSM objects, a dictionary of GPL objects and its own metadata::

    # imports
    import GEOparse
    import pandas as pd

    # prepare data for columns and table
    columns_header = ['a', 'b', 'c']
    columns_description = ['first column', 'second column', 'third column']
    table_data = [[1, 2, 3],
                 [4, 5, 6]]

    # create table and columns
    table = pd.DataFrame(table_data, columns=columns_header)
    columns = pd.DataFrame(columns_description, columns_header)
    columns.columns = ['description'] # columns header should contain description

    # prepare metadata for objects. Each value of the dictionary should be a list
    gpl_metadata = {'name': ['FooGPL']}
    gsm_metadata = {'name': ['FooGSM']}
    metadata = {'name': ['FooGSE']}

    # initialize GPL and GSM object(s)
    gpl = GEOparse.GPL(name='FooGPL', table=table, metadata=gpl_metadata, columns=columns)
    gsm = GEOparse.GSM(name='FooGSM', table=table, metadata=gsm_metadata, columns=columns)

    # prepare attributes for GSE
    gsms = {'FooGSM': gsm}
    gpls = {'FooGPL': gpl}

    # initialize GSE
    gse = GEOparse.GSE(name='FooGSE', metadata=metadata, gpls=gpls, gsms=gsms)

    # save gse as SOFT file
    gse.to_soft("./GSEFoo.soft")

This creates file with following content::

    ^SERIES = FooGSE
    !Series_name = FooGSE
    ^SAMPLE = FooGSM
    !Sample_name = FooGSM
    #a = first column
    #b = second column
    #c = third column
    !sample_table_begin
    a   b   c
    1   2   3
    4   5   6
    !sample_table_end
    ^PLATFORM = FooGPL
    !Platform_name = FooGPL
    #a = first column
    #b = second column
    #c = third column
    !platform_table_begin
    a   b   c
    1   2   3
    4   5   6
    !platform_table_end

Of course in this case the file is simpler than the code that generates it but in normal situation
this is reversed. For more information on what GEO objects are available and what parameters one need
to create them see the :ref:`working-with-objects` section.

.. _working-with-objects:
Working with GEO objects
^^^^^^^^^^^^^^^^^^^^^^^^

BaseGEO
-------

All GEO objects inherit from abstract base class :class:`GEOparse.BaseGEO`. Two main attributes of that
class are the *name* and *metadata*. 

*metadata* is a dictionary of useful information about samples which occurs in the SOFT file with bang (!)
in the beginning. Each value of this dictionary id a list (even with one element).

GSM (Sample)
------------

A GSM (or a Sample) contains information the conditions and preparation of a Sample. In the GEO database
sample is assigned to unique and stable GEO accession number that is composed of  'GSM' followed by numbers eg. GSM906.

In GEOparse Sample is represented by GEOparse.GSM object that contains tree main attributes:
 * inherited from BaseGEO :attr:`metadata`
 * :attr:`table` -- :class:`pandas.DataFrame` with the data table from SOFT file
 * :attr:`columns` -- :class:`pandas.DataFrame` that contains *description* column with the
   information about columns in the :attr:`table`

See API for more information.

GPL (Platform)
--------------

A GPL (or a Platform) contains a tab-delimited table containing the array definition eg. mappings from probe IDs to
RefSeq IDs. Similarly to GSM, it is assigned to unique and stable GEO accession number that is composed of  'GPL'
followed by numbers eg. GPL2020.

In GEOparse Platform is represented by GEOparse.GSM object that contains tree main attributes:
 * inherited from BaseGEO :attr:`metadata`
 * :attr:`table` -- :class:`pandas.DataFrame` with the data table from SOFT file
 * :attr:`columns` -- :class:`pandas.DataFrame` that contains *description* column with the
   information about columns in the :attr:`table`

See API for more information.

GSE (Series)
------------

A GSE (or a Series) is an original submitter-supplied record that summarizes whole study including samples and platforms.
GSE is assigned to unique and stable GEO accession number that starts at GSE followed by numbers eg. GSE1563.

In GEOparse Series is represented by GEOparse.GSE object that contains tree main attributes:
 * inherited from BaseGEO :attr:`metadata`
 * :attr:`gsms` -- :class:`dict` with all samples in this GSE as GSM objects
 * :attr:`gpls` -- :class:`dict` with all platforms in this GSE as GSM objects

See API for more information.

GDS (Dataset)
-------------

A GDS (or a Dataset) is a curated file that hold a summarised combination of a Series file and its samples.
GDS is assigned to unique and stable GEO accession number that starts at GDS followed by numbers eg. GDS1563.

In GEOparse Dataset is represented by GEOparse.GDS object that contains tree main attributes:
 * inherited from BaseGEO :attr:`metadata`
 * :attr:`table` -- :class:`pandas.DataFrame` with the data table from SOFT file
 * :attr:`columns` -- :class:`pandas.DataFrame` that contains *description* column with the
   information about columns in the :attr:`table` and additional information according to GDS file

See API for more information.

Examples
^^^^^^^^

.. toctree::
    :maxdepth: 2
    :titlesonly:
    :glob:
    :hidden:

    Analyse_hsa-miR-124a-3p_transfection_time-course.rst
