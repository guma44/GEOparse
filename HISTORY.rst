.. :changelog:

History
-------

0.1.8 (2016-11-02)
---------------------

Thanks to Tycho Bismeijer:

* Python 3 Compatibility
* Bio.Entrez dependency optional


0.1.7 (2016-05-30)
---------------------

Thanks to Simon van Heeringen:


* bugfix in datasets with multiple associated relations
* --split-files to fastq-dump to support paired-end experiments by default
* parse a GPL that also contains series and sample information
* gsm2fastq command to make download easier
* initial Aspera download support


0.1.6 (2016-04-12)
---------------------

* Bugfixes
* SRA function of GSE can now filter GSMs


0.1.5 (2016-02-03)
---------------------

* Added functions to download supplementary files including raw files from SRA

0.1.4 (2015-09-27)
---------------------

* Updated documentation including example
* Updated tests: they now cover 80% of library with all important functions
* Added pivot_and_annotate method to GSE object
* Bugfixes

0.1.3 (2015-08-30)
---------------------

* Updated documentation
* Added pivot_samples to GSE object
* Code of GEOTypes was refactored
* All objects now have to_soft function
* Various bugfixes

0.1.2 (2015-08-23)
---------------------

* Added GDS support
* Added to_soft methods to GSE, GSM and GPL
* Added DATABASE entry support to GSE and GDS

0.1.1 (2015-08-16)
---------------------

* Brown-Bag release

0.1.0 (2015-08-16)
---------------------

* First release on PyPI.
