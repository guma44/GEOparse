import os
import sys
import gzip
import glob
from errno import EEXIST
from contextlib import closing
from shutil import copyfileobj
from contextlib import contextmanager

try:
    from urllib.request import urlopen
    from urllib.error import URLError
except ImportError:
    from urllib2 import urlopen, URLError
import subprocess as sp
from six import iteritems

from .downloader import Downloader
from .logger import geoparse_logger as logger


def mkdir_p(path_to_dir):
    """Make directory(ies).

    This function behaves like mkdir -p.

    Args:
        path_to_dir (:obj:`str`): Path to the directory to make.
    """
    try:
        os.makedirs(path_to_dir)
    except OSError as e:  # Python >2.5
        if e.errno == EEXIST and os.path.isdir(path_to_dir):
            logger.debug(
                "Directory %s already exists. Skipping." % path_to_dir)
        else:
            raise e


def download_from_url(url, destination_path,
                      force=False, aspera=False, silent=False):
    """Download file from remote server.

    If the file is already downloaded and  ``force`` flag is on the file will
    be removed.

    Args:
        url (:obj:`str`): Path to the file on remote server (including file
            name)
        destination_path (:obj:`str`): Path to the file on local machine
            (including file name)
        force (:obj:`bool`): If file exist force to overwrite it. Defaults to
            False.
        aspera (:obj:`bool`): Download with Aspera Connect. Defaults to False.
        silent (:obj:`bool`): Do not print any message. Defaults to False.
    """
    if aspera and url.startswith("http"):
        logger.warning("Aspera Connect allows only FTP servers - falling "
                       "back to normal download")
        aspera = False
    use_http_for_ftp = os.environ.get("GEOPARSE_USE_HTTP_FOR_FTP") == "yes"
    if os.environ.get("http_proxy") is not None or use_http_for_ftp:
        if url.startswith("ftp://"):
            url = url.replace("ftp://", "http://")
            logger.warning("Changing FTP to HTTP: %s" % url)
    try:
        fn = Downloader(
            url,
            outdir=os.path.dirname(destination_path),
            filename=os.path.basename(destination_path))
        if aspera:
            fn.download_aspera(
                user="anonftp",
                host="ftp-trace.ncbi.nlm.nih.gov",
                silent=silent)
        else:
            fn.download(silent=silent, force=force)
    except URLError:
        logger.error("Cannot find file %s" % url)


@contextmanager
def smart_open(filepath):
    """Open file intelligently depending on the source and python version.

    Args:
        filepath (:obj:`str`): Path to the file.

    Yields:
        Context manager for file handle.

    """
    if filepath[-2:] == "gz":
        mode = "rt"
        fopen = gzip.open
    else:
        mode = "r"
        fopen = open
    if sys.version_info[0] < 3:
        fh = fopen(filepath, mode)
    else:
        fh = fopen(filepath, mode, errors="ignore")
    try:
        yield fh
    except IOError:
        fh.close()
    finally:
        fh.close()


def which(program):
    """Check if executable exists.

    The code is taken from:
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    Args:
        program (:obj:`str`): Path to the executable.

    Returns:
        :obj:`str` or :obj:`None`: Path to the program or None.

    """

    def is_exe(fpath):
        """Check if fpath is executable."""
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
