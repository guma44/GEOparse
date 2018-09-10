===============================
GEOparse
===============================

.. image:: https://img.shields.io/pypi/v/GEOparse.svg
        :target: https://pypi.python.org/pypi/GEOparse

.. image:: https://img.shields.io/travis/guma44/GEOparse.svg
        :target: https://travis-ci.org/guma44/GEOparse


Python library to access Gene Expression Omnibus Database (GEO).

GEOparse is python package that can be used to query and retrieve data from Gene Expression Omnibus database (GEO).
The inspiration and the base for it is great R library GEOquery.

* Free software: BSD license
* Documentation: https://GEOparse.readthedocs.org.

Features
--------

* Download GEO series, datasets etc. as SOFT files
* Download supplementary files for GEO series to use them locally
* Load GEO SOFT as easy to use and manipulate objects
* Prepare your data for GEO upload

Installation
------------

At the command line::

    $ pip install GEOparse

TODO
----

There is still work to do so any contribution is welcome. Any bug/error that you report
will improve the library.

The main issues are:

* add checking for compatibility with SOFT files
* expand GEOTypes objects with useful functions for differential expression analysis
* share your idea
* add more tests - that's always good idea :)
