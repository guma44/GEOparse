.. :changelog:

History
-------

2.0.3 (2021-02-03)
------------------

 * BUGFIX: shutil.move can lead to "operation not permitted" error
 * Allow to download small GSE object by using `how` parameter



2.0.2 (2020-08-26)
------------------

 * Allow to pass kwargs to open function to controll eg. encoding

2.0.1 (2020-02-20)
------------------

 * prevent installation on old python versions

2.0.0 (2020-02-20)
------------------

 * Drop support for python 2

 * Fix issue with SRR entiries not being available through FTP server

 * Better error message if entry is not available

 * Better code formatting

 * Some small improvements


1.2.0 (2019-08-21)
------------------

 * Remove BioPython dependency - use requests instead

 * Replace wget with custom file downloader

 * Allow to parse GSE file partially

 * Support for parallel-fastq-dump

 * Allow to use proxy server (FTP)

 * Fix Travis-CI issues with FTP download

 * Various bugfixes

1.1.0 (2018-09-10)
------------------

 * Allow parallel download of supplementary files
 
 * Add timeout to urlopen callsq

 * Bugfix: Keep `logger` as a module name and rename `logger` object to
   `geoparse_logger`. This is a breaking change.

 * Bugfix: Some samples do not have table data error in python 3
 
 * Bugfix: broken download when supplementary_files is empty or contains invalid URLs
 
 * Some small bugfixes
 

1.0.5 (2018-01-12)
------------------

 * Bugfix: Some samples do not have table data error in python 3

1.0.4 (2018-01-08)
------------------

 * Bugfix: Empty line in the SOFT file caused an error in python 3

1.0.3 (2017-11-01)
------------------

 * Bugfix: Fixed the FTP link

1.0.2 (2017-11-01)
------------------

 * Bugfix: type name was depended on the order of entries

1.0.1 (2017-08-10)
------------------

 * Hotfix: wrong path split in Windows

1.0.0 (2017-07-21)
------------------

* Many small bug fixes:
  * unknown subset types added to columns
  * silent=True is really silent
  * correct treatment of duplicated columns
  * illegal file names and no filtering of user input from GEO to create the file names
  * platform was not imported but used
  * fixed issues of python 2 and 3 compatibility
* Logging replaced stdout and stderr + ability to set verbosity and log file
* Return downloaded paths from download functions
* Updated documentation according to Google docstring style guide
* Tests update
* Code refactored to be more PEP-8 friendly


0.1.10 (2017-03-27)
-------------------

* Important fix for SRA download
* Fix duplicated columns issue
* Python 2 and 3 open compatibility


0.1.9 (2017-03-10)
------------------

* Added property phenotype_data to access phenotype data of GSE
* Fixed windows issue with file names
* replaced default download function with wgetter
* Update documentation
* Various bugfixes

0.1.8 (2016-11-02)
------------------

Thanks to Tycho Bismeijer:

* Python 3 Compatibility
* Bio.Entrez dependency optional


0.1.7 (2016-05-30)
------------------

Thanks to Simon van Heeringen:


* bugfix in datasets with multiple associated relations
* --split-files to fastq-dump to support paired-end experiments by default
* parse a GPL that also contains series and sample information
* gsm2fastq command to make download easier
* initial Aspera download support


0.1.6 (2016-04-12)
------------------

* Bugfixes
* SRA function of GSE can now filter GSMs


0.1.5 (2016-02-03)
------------------

* Added functions to download supplementary files including raw files from SRA

0.1.4 (2015-09-27)
------------------

* Updated documentation including example
* Updated tests: they now cover 80% of library with all important functions
* Added pivot_and_annotate method to GSE object
* Bugfixes

0.1.3 (2015-08-30)
------------------

* Updated documentation
* Added pivot_samples to GSE object
* Code of GEOTypes was refactored
* All objects now have to_soft function
* Various bugfixes

0.1.2 (2015-08-23)
------------------

* Added GDS support
* Added to_soft methods to GSE, GSM and GPL
* Added DATABASE entry support to GSE and GDS

0.1.1 (2015-08-16)
------------------

* Brown-Bag release

0.1.0 (2015-08-16)
------------------

* First release on PyPI.
