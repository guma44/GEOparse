import glob
import gzip
import os
import subprocess as sp
import sys
import time
from contextlib import closing, contextmanager
from errno import EEXIST
from random import choice
from shutil import copyfileobj

import requests
from six import iteritems

from .downloader import Downloader
from .logger import geoparse_logger as logger

try:
    from urllib.request import urlopen
    from urllib.error import URLError
except ImportError:
    from urllib2 import urlopen, URLError


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
            logger.debug("Directory %s already exists. Skipping." % path_to_dir)
        else:
            raise e


def download_from_url(url, destination_path, force=False, aspera=False, silent=False):
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
        logger.warning(
            "Aspera Connect allows only FTP servers - falling "
            "back to normal download"
        )
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
            filename=os.path.basename(destination_path),
        )
        if aspera:
            fn.download_aspera(
                user="anonftp", host="ftp-trace.ncbi.nlm.nih.gov", silent=silent
            )
        else:
            fn.download(silent=silent, force=force)
    except Exception as err:
        raise IOError(
            "Download failed due to '%s'. ID could be incorrect or the " % err
            + "data might not be public yet."
        )


@contextmanager
def smart_open(filepath, **open_kwargs):
    """Open file intelligently depending on the source and python version.

    Args:
        filepath (:obj:`str`): Path to the file.

    Yields:
        Context manager for file handle.

    """
    if "errors" not in open_kwargs:
        open_kwargs["errors"] = "ignore"
    if filepath[-2:] == "gz":
        open_kwargs["mode"] = "rt"
        fopen = gzip.open
    else:
        open_kwargs["mode"] = "r"
        fopen = open
    if sys.version_info[0] < 3:
        fh = fopen(filepath, **open_kwargs)
    else:
        fh = fopen(filepath, **open_kwargs)
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


def esearch(db, term, **kwargs):
    """Run an Entrez search.

    ESearch searches and retrieves primary IDs (for use in EFetch, ELink
    and ESummary) and term translations.
    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch

    Arguments:
        db {:class:`str`} -- A DB to use
        term {:class:`str`} -- term to search

    Returns:
        :class:`dict` -- results of the search
    """
    time.sleep(choice(range(10)))
    wait_time = 10
    number_of_trials = 10
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    data = {
        "db": db,
        "term": term,
        "usehistory": "y",
        "retmode": "json",
        "tool": "geoparse",
    }
    data.update(kwargs)
    is_successfull = False
    error = None
    for trial in range(number_of_trials):
        try:
            res = requests.post(url, data=data)
            res.raise_for_status()
            is_successfull = True
            break
        except requests.exceptions.HTTPError as httperr:
            error = httperr
            logger.warning(
                "An error occurred: %s for %s with data %s" % (httperr, url, str(data))
            )
            if httperr.response.status_code == 429:
                # This means that there is too many requests
                try:
                    header_wait_time = int(httperr.headers["Retry-After"])
                except Exception:
                    header_wait_time = wait_time
                logger.warning(
                    ("%s, trial %i out of %i, waiting " "for %i seconds.")
                    % (str(httperr), trial, number_of_trials, header_wait_time)
                )
                time.sleep(header_wait_time)
            else:
                logger.warning(
                    ("%s, trial %i out of %i, waiting " "for %i seconds.")
                    % (str(httperr), trial, number_of_trials, wait_time)
                )
                time.sleep(wait_time)
    if not is_successfull:
        if error is not None:
            raise error
        else:
            raise RuntimeError("Failed in esearch for unknown reason")
    return res.json()


def efetch(db, **kwargs):
    """Fetch Entrez results.

    Retrieves records in the requested format from a list of one or
    more UIs.
    See the online documentation for an explanation of the parameters:
    http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

    Args:
        db {:class:`str`} -- A DB to use. eg. sra

    Returns:
        :class:`dict`-- results of the fetch
    """
    time.sleep(choice(range(10)))
    number_of_trials = 10
    wait_time = 10
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    data = {"db": db, "retmode": "text", "tool": "geoparse"}
    data.update(kwargs)
    try:
        ids = data["id"]
    except KeyError:
        pass
    else:
        if isinstance(ids, list):
            ids = ",".join(ids)
            data["id"] = ids
        elif isinstance(ids, int):
            ids = str(ids)
            data["id"] = ids
    is_successfull = False
    error = None
    for trial in range(number_of_trials):
        try:
            res = requests.post(url, data=data)
            res.raise_for_status()
            is_successfull = True
            break
        except requests.exceptions.HTTPError as httperr:
            error = httperr
            logger.warning(
                "An error occurred: %s for %s with data %s" % (httperr, url, str(data))
            )
            if httperr.response.status_code == 429:
                # This means that there is too many requests
                try:
                    header_wait_time = int(httperr.headers["Retry-After"])
                except Exception:
                    header_wait_time = wait_time
                logger.warning(
                    ("%s, trial %i out of %i, waiting " "for %i seconds.")
                    % (str(httperr), trial, number_of_trials, header_wait_time)
                )
                time.sleep(header_wait_time)
            else:
                logger.warning(
                    ("%s, trial %i out of %i, waiting " "for %i seconds.")
                    % (str(httperr), trial, number_of_trials, wait_time)
                )
                time.sleep(wait_time)
    if not is_successfull:
        if error is not None:
            raise error
        else:
            raise RuntimeError("Failed in efetch for unknown reason")

    return res.text
