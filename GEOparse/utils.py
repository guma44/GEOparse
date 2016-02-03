import os
from errno import EEXIST
from sys import stderr, stdout
from contextlib import closing
from shutil import copyfileobj
from urllib2 import urlopen, URLError

def mkdir_p(path_to_dir):
    try:
        os.makedirs(path_to_dir)
    except OSError as e: # Python >2.5
        if e.errno == EEXIST and os.path.isdir(path_to_dir):
            stderr.write("Directory %s already exists. Skipping.\n" % path_to_dir)
        else:
            raise e


def download_from_url(url, destination_path, force=False):
    """Download file from remote server

    :param url: path to the file on remote server (including file name)
    :param destination_path: path to the file on local machine (including file name)
    :param force: bool - if file exist force to overwrite it , defaults to False
    """
    try:
        is_already_downloaded = os.path.isfile(destination_path)
        if is_already_downloaded:
            if force:
                with closing(urlopen(url)) as r:
                    with open(destination_path, mode='wb') as f:
                        stderr.write("Downloading %s to %s\n" % (url, destination_path))
                        copyfileobj(r, f)
            else:
                stderr.write("File already exist. Use force=True if you would like to overwrite it.\n")
        else:
            with closing(urlopen(url)) as r:
                with open(destination_path, mode='wb') as f:
                    stderr.write("Downloading %s to %s\n" % (url, destination_path))
                    copyfileobj(r, f)
    except URLError:
        stderr.write("Cannot find file %s" % url)
