import logging
from logging.handlers import RotatingFileHandler

########################################################################
# Set up logging
########################################################################

_formatter = logging.Formatter(
    fmt='%(asctime)s %(levelname)s %(module)s - %(message)s',
    datefmt="%d-%b-%Y %H:%M:%S")
_console_handler = logging.StreamHandler()
_console_handler.setFormatter(_formatter)
logger = logging.getLogger('GEOparse')
logger.setLevel(logging.getLevelName('DEBUG'))
logger.addHandler(_console_handler)


def set_verbosity(level):
    """Set the log level.

    Args:
        level (:obj:`str`): Level name eg. DEBUG or ERROR
    """
    logger.setLevel(logging.getLevelName(level))


def add_log_file(path):
    """Add log file.

    Args:
        path (:obj:`str`): Path to the log file.
    """
    logfile_handler = RotatingFileHandler(
        path, maxBytes=50000, backupCount=2)
    formatter = logging.Formatter(
        fmt='%(asctime)s %(levelname)s %(module)s - %(message)s',
        datefmt="%d-%b-%Y %H:%M:%S")
    logfile_handler.setFormatter(formatter)
    logger.addHandler(logfile_handler)


__all__ = ('logger', 'set_verbosity', 'add_log_file')
