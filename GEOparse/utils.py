import os
import sys
from errno import EEXIST
from sys import stderr, stdout
from contextlib import closing
from shutil import copyfileobj
try:
    from urllib.request import urlopen
    from urllib.error import URLError
except ImportError:
    from urllib2 import urlopen, URLError
import subprocess as sp
import wgetter

def mkdir_p(path_to_dir):
    try:
        os.makedirs(path_to_dir)
    except OSError as e: # Python >2.5
        if e.errno == EEXIST and os.path.isdir(path_to_dir):
            stderr.write("Directory %s already exists. Skipping.\n" % path_to_dir)
        else:
            raise e

def download_aspera(url, dest_path, user="anonftp", ftp="ftp-trace.ncbi.nlm.nih.gov"):
    sys.stderr.write("Downloading {} using aspera\n".format(url))
    aspera_home = os.environ.get("ASPERA_HOME", None)
    if not aspera_home:
        raise ValueError("environment variable $ASPERA_HOME not set")
    if not os.path.exists(aspera_home):
        raise ValueError("$ASPERA_HOME directory {} does not exist".format(aspera_home))
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
    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()

def download_from_url(url, destination_path, force=False, aspera=False):
    """Download file from remote server

    :param url: path to the file on remote server (including file name)
    :param destination_path: path to the file on local machine (including file name)
    :param force: bool - if file exist force to overwrite it , defaults to False
    """
    try:
        is_already_downloaded = os.path.isfile(destination_path)
        if is_already_downloaded:
            if force:
                stderr.write("Downloading %s to %s\n" % (url, destination_path))
                fn = wgetter.download(url, outdir=os.path.dirname(destination_path))
            else:
                stderr.write("File already exist. Use force=True if you would like to overwrite it.\n")
        else:
            if aspera:
                download_aspera(url, destination_path)
            else:
                stderr.write("Downloading %s to %s\n" % (url, destination_path))
                fn = wgetter.download(url, outdir=os.path.dirname(destination_path))
    except URLError:
        stderr.write("Cannot find file %s" % url)
