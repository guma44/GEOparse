import os
import requests

from tqdm import tqdm
from ftplib import FTP
from urlparse import urlparse

from .logger import geoparse_logger as logger

def get_filename(url):
    filename = os.path.basename(urlparse(url).path).strip(" \n\t.")
    if len(filename) == 0:
        raise Exception("Cannot parse filename from %s" % url)
    return filename

def download(url, outdir, filename=None, silent=False):
    """Download from URL.

    :param url: URL - could be FTP or HTTP
    :type url: :class:`str`
    :param outdir: Path to the outdir file
    :type outdir: :class:`str`
    :raises: ValueError if URL is invalid
    """
    if filename is None:
        filename = get_filename(url)
    destination = os.path.join(outdir, filename)
    logger.info("Downloading %s to %s" % (url, destination))
    if url.startswith("http"):
        r = requests.get(url, stream=True)
        # Total size in bytes.
        total_size = int(r.headers.get('content-length', 0))
        chunk_size = 1024
        with open(destination, 'wb') as f:
            if silent:
                for chunk in r.iter_content(chunk_size):
                    if chunk:
                        f.write(chunk)
            else:
                with tqdm(
                        total=total_size,
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                        leave=True) as pbar:
                    for chunk in r.iter_content(chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
    elif url.startswith("ftp"):
        parsed_url = urlparse(url)
        try:
            ftp = FTP(parsed_url.netloc)
            ftp.login()
            total_size = ftp.size(parsed_url.path)
            if total_size is None:
                total_size = 0
            with open(destination, 'wb') as f:
                with tqdm(total=total_size,
                          unit="B",
                          unit_scale=True,
                          unit_divisor=1024,
                          leave=True) as pbar:
                    def _write(data):
                        pbar.update(len(data))
                        f.write(data)
                    ftp.retrbinary("RETR %s" % parsed_url.path, _write)
            ftp.quit()
        except Exception:
            try:
                ftp.quit()
                logger.error("Error when trying to retreive %s." % url,
                             exc_info=True)
            except Exception:
                logger.error("Error when quiting FTP server.", exc_info=True)

    else:
        raise ValueError("Invalid URL %s" % url)
