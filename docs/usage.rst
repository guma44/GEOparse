========
Usage
========

To simplest example of usage::

    import GEOparse

    gse = GEOparse.get_GEO(geo="GSE1563", destdir="./")

    print
    print "GSM example:"
    for gsm_name, gsm in gse.gsms.iteritems():
        print "Name: ", gsm_name
        print "Metadata:",
        for key, value in gsm.metadata.iteritems():
            print " - %s : %s" % (key, ", ".join(value))
        print "Table data:",
        print gsm.table.head()
        break

    print
    print "GPL example:"
    for gpl_name, gpl in gse.gpls.iteritems():
        print "Name: ", gpl_name
        print "Metadata:",
        for key, value in gpl.metadata.iteritems():
            print " - %s : %s" % (key, ", ".join(value))
        print "Table data:",
        print gpl.table.head()
        break

Working with GEO accession
^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Working with SOFT files
^^^^^^^^^^^^^^^^^^^^^^^

TODO

Working with GEO objects
^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Examples
^^^^^^^^

TODO
