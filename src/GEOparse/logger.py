import logging
from logging.handlers import RotatingFileHandler

########################################################################
# Set up logging
########################################################################

_formatter = logging.Formatter(
    fmt="%(asctime)s %(levelname)s %(module)s - %(message)s",
    datefmt="%d-%b-%Y %H:%M:%S",
)
_console_handler = logging.StreamHandler()
_console_handler.setFormatter(_formatter)
geoparse_logger = logging.getLogger("GEOparse")
geoparse_logger.setLevel(logging.getLevelName("DEBUG"))
geoparse_logger.addHandler(_console_handler)


def set_verbosity(level):
    """Set the log level.

    Args:
        level (:obj:`str`): Level name eg. DEBUG or ERROR
    """
    geoparse_logger.setLevel(logging.getLevelName(level))


def add_log_file(path):
    """Add log file.

    Args:
        path (:obj:`str`): Path to the log file.
    """
    logfile_handler = RotatingFileHandler(path, maxBytes=50000, backupCount=2)
    formatter = logging.Formatter(
        fmt="%(asctime)s %(levelname)s %(module)s - %(message)s",
        datefmt="%d-%b-%Y %H:%M:%S",
    )
    logfile_handler.setFormatter(formatter)
    geoparse_logger.addHandler(logfile_handler)


__all__ = ("geoparse_logger", "set_verbosity", "add_log_file")
