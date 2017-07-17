# -*- coding: utf-8 -*-

__author__ = 'Rafal Gumienny'
__email__ = 'guma44@gmail.com'
__version__ = '0.1.10'

from .GEOparse import get_GEO, get_GEO_file, parse_GPL, parse_GSE, parse_GSM
from .GEOTypes import (DataIncompatibilityException,
                       NoMetadataException,
                       GSM,
                       GPL,
                       GDSSubset,
                       GEODatabase,
                       GDS,
                       GSE)

import logging
from logging.handlers import RotatingFileHandler
from .logger import logger


def setVerbosity(level):
    """Set the log level."""
    logger.setLevel(logging.getLevelName(level))


def addLogFile(path):
    logfile_handler = RotatingFileHandler(
        path, maxBytes=50000, backupCount=2)
    formatter = logging.Formatter(
        fmt='%(asctime)s %(levelname)s %(module)s - %(message)s',
        datefmt="%d-%b-%Y %H:%M:%S")
    logfile_handler.setFormatter(formatter)
    logger.addHandler(logfile_handler)
