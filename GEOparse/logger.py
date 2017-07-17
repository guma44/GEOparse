import logging


########################################################################
# Set up logging
########################################################################

_formatter = logging.Formatter(fmt='%(asctime)s %(levelname)s %(module)s - %(message)s',
                               datefmt="%d-%b-%Y %H:%M:%S")
_console_handler = logging.StreamHandler()
_console_handler.setFormatter(_formatter)
logger = logging.getLogger('GEOparse')
logger.setLevel(logging.getLevelName('DEBUG'))
logger.addHandler(_console_handler)


__all__ = ('logger',)
