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
import wgetter
from six import iteritems

from .logger import logger


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


def download_aspera(url, dest_path,
                    user="anonftp",
                    ftp="ftp-trace.ncbi.nlm.nih.gov"):
    """Download file with Aspera Connect.

    For details see the documentation ov Aspera Connect

    Args:
        url (:obj:`str`): URL to the file
        dest_path (:obj:`str`): Destination path.
        user (:obj:`str`, optional): User. Defaults to anonftp.
        ftp (:obj:`str`, optional): FTP path. Defaults to
            "ftp-trace.ncbi.nlm.nih.gov".
    """
    logger.info("Downloading {} using aspera\n".format(url))
    aspera_home = os.environ.get("ASPERA_HOME", None)
    if not aspera_home:
        raise ValueError("environment variable $ASPERA_HOME not set")
    if not os.path.exists(aspera_home):
        raise ValueError(
            "$ASPERA_HOME directory {} does not exist".format(aspera_home))
    ascp = os.path.join(aspera_home, "connect/bin/ascp")
    key = os.path.join(aspera_home, "connect/etc/asperaweb_id_dsa.openssh")
    if not os.path.exists(ascp):
        raise ValueError("could not find ascp binary")
    if not os.path.exists(key):
        raise ValueError("could not find openssh key")

    if url.startswith("ftp://"):
        url = url.replace("ftp://", "")
    url = url.replace(ftp, "")

    cmd = "{} -i {} -k1 -T -l400m {}@{}:{} {}".format(
        ascp, key, user, ftp, url, dest_path)
    logger.debug(cmd)
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()


def download_from_url(url, destination_path,
                      force=False, aspera=False, silent=False):
    """Download file from remote server

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
    try:
        is_already_downloaded = os.path.isfile(destination_path)
        if is_already_downloaded:
            if force:
                if not silent:
                    logger.info(
                        "Downloading %s to %s" % (url, destination_path))
                    fn = wgetter.download(url, outdir=os.path.dirname(
                        destination_path))
                else:
                    with closing(urlopen(url)) as r:
                        with open(destination_path, mode='wb') as f:
                            copyfileobj(r, f)
            else:
                logger.info("File already exist. Use force=True if you would "
                            "like to overwrite it.")
        else:
            if aspera:
                download_aspera(url, destination_path)
            else:
                if not silent:
                    logger.info(
                        "Downloading %s to %s" % (url, destination_path))
                    fn = wgetter.download(url, outdir=os.path.dirname(
                        destination_path))
                else:
                    with closing(urlopen(url)) as r:
                        with open(destination_path, mode='wb') as f:
                            copyfileobj(r, f)
    except URLError:
        logger.error("Cannot find file %s" % url)


def download_unpack_SRA_for_parallel(args):
    """
    Auxiliary function for parallel download of sra.
    """
    return download_unpack_SRA(*args)


def download_unpack_SRA(path, ftpaddres, directory_path,
                        filetype='fasta',
                        force=False,
                        aspera=False,
                        silent=False,
                        fastq_dump_options=None,
                        keep_sra=False):
    """
    Combination of download_from_url for sra and unpacking with fastq-dump.

    :param path: downloaded path
    :param ftpaddres: ftp address
    :param directory_path: target local directory
    :param filetype: 'fastq' or 'fasta' for fastq-dump, 'sra' to leave SRA unpacked
    :param force: overwrite existing sra?
    :param aspera: download with aspera
    :param silent: supress wgetter log (get rid of enormous log file)
    :param fastq_dump_options: options for fastq-dump, see .download_SRA description
    :param keep_sra: keep original sra for later use
    :return: downloaded paths (note that if sequencing is pair-ended it might generate list of output files)
    """

    # Make the directory

    mkdir_p(os.path.abspath(directory_path))

    sra_run = path.split("/")[-1]
    logger.info("Analysing %s" % sra_run)
    url = ftpaddres.format(range_subdir=sra_run[:6],
                           file_dir=sra_run)
    logger.debug("URL: %s", url)
    filepath = os.path.abspath(
        os.path.join(directory_path, "%s.sra" % sra_run))
    download_from_url(url, filepath, aspera=aspera, silent=silent, force=force)

    if filetype in ["fasta", "fastq"]:
        if which('fastq-dump') is None:
            logger.error("fastq-dump command not found")
        ftype = ""
        if filetype == "fasta":
            ftype = " --fasta "
        cmd = "fastq-dump"
        for fqoption, fqvalue in iteritems(fastq_dump_options):
            if fqvalue:
                cmd += (" --%s %s" % (fqoption, fqvalue))
            else:
                cmd += (" --%s" % fqoption)
        cmd += " %s --outdir %s %s"
        cmd = cmd % (ftype, directory_path, filepath)
        logger.debug(cmd)
        process = sp.Popen(cmd, stdout=sp.PIPE,
                           stderr=sp.PIPE,
                           shell=True)
        logger.info("Converting to %s/%s*.%s.gz\n" % (
            directory_path, sra_run, filetype))
        pout, perr = process.communicate()
        downloaded_path = glob.glob(os.path.join(
            directory_path,
            "%s*.%s.gz" % (sra_run, filetype)
        ))

    elif filetype == 'sra':
        downloaded_path = glob.glob(os.path.join(
            directory_path,
            "%s*.%s" % (sra_run, filetype)
        ))

    else:
        downloaded_path = glob.glob(os.path.join(
            directory_path,
            "%s*" % sra_run))
        logger.error("Filetype %s not supported." % filetype)

    if not keep_sra and filetype != 'sra':
        # Delete sra file
        os.unlink(filepath)

    return downloaded_path


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
