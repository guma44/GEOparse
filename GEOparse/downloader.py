import os
import requests

from tqdm import tqdm
from ftplib import FTP
try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

from .logger import geoparse_logger as logger

class Downloader(object):
    """Downloader class."""
    def __init__(self, url, outdir, filename=None, silent=False):
        self.url = url
        if outdir is None:
            self.outdir = os.getcwd()
        else:
            self.outdir = outdir
        if filename is None:
            self.filename = self._get_filename()
        else:
            self.filename = filename
        self.silent = silent

    @property
    def destination(self):
        """Get the destination path.

        This is the property should be calculated every time it is used because
        a user could change the outdir and filename dynamically.
        """
        return os.path.join(os.path.abspath(self.outdir), self.filename)

    def download(self):
        """Download from URL."""
        logger.info("Downloading %s to %s" % (self.url, self.destination))
        if self.url.startswith("http"):
            self._download_http()
        elif self.url.startswith("ftp"):
            self._download_ftp()
        else:
            raise ValueError("Invalid URL %s" % self.url)

    def _get_filename(self):
        filename = os.path.basename(urlparse(self.url).path).strip(" \n\t.")
        if len(filename) == 0:
            raise Exception("Cannot parse filename from %s" % self.url)
        return filename

    def _download_ftp(self):
        parsed_url = urlparse(self.url)
        try:
            ftp = FTP(parsed_url.netloc)
            ftp.login()
            total_size = ftp.size(parsed_url.path)
            if total_size is None:
                total_size = 0
            with open(self.destination, 'wb') as f:
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
                logger.error("Error when trying to retreive %s." % self.url,
                             exc_info=True)
            except Exception:
                logger.error("Error when quiting FTP server.", exc_info=True)

    def _download_http(self):
        r = requests.get(self.url, stream=True)
        # Total size in bytes.
        total_size = int(r.headers.get('content-length', 0))
        chunk_size = 1024
        with open(self.destination, 'wb') as f:
            if self.silent:
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
