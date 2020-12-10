"""A module used for downloading files."""
import hashlib
import os
import shutil
import subprocess as sp
import tempfile
from ftplib import FTP

import requests
from tqdm import tqdm

from .logger import geoparse_logger as logger

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse


class Downloader(object):
    """Downloader class."""

    def __init__(self, url, outdir=None, filename=None):
        self.url = url
        if outdir is None:
            self.outdir = os.getcwd()
        else:
            self.outdir = outdir
        if filename is None:
            self.filename = self._get_filename()
        else:
            self.filename = filename

        with tempfile.NamedTemporaryFile(delete=True) as tmpf:
            self._temp_file_name = tmpf.name

    @property
    def destination(self):
        """Get the destination path.

        This is the property should be calculated every time it is used because
        a user could change the outdir and filename dynamically.
        """
        return os.path.join(os.path.abspath(self.outdir), self.filename)

    def download(self, force=False, silent=False):
        """Download from URL."""

        def _download():
            if self.url.startswith("http"):
                self._download_http(silent=silent)
            elif self.url.startswith("ftp"):
                self._download_ftp(silent=silent)
            else:
                raise ValueError("Invalid URL %s" % self.url)
            logger.debug("Moving %s to %s" % (self._temp_file_name, self.destination))
            shutil.copyfile(self._temp_file_name, self.destination)
            logger.debug("Successfully downloaded %s" % self.url)

        try:
            is_already_downloaded = os.path.isfile(self.destination)
            if is_already_downloaded:
                if force:
                    try:
                        os.remove(self.destination)
                    except Exception:
                        logger.error("Cannot delete %s" % self.destination)
                    logger.info("Downloading %s to %s" % (self.url, self.destination))
                    logger.debug(
                        "Downloading %s to %s" % (self.url, self._temp_file_name)
                    )
                    _download()
                else:
                    logger.info(
                        (
                            "File %s already exist. Use force=True if you"
                            " would like to overwrite it."
                        )
                        % self.destination
                    )
            else:
                _download()
        finally:
            try:
                os.remove(self._temp_file_name)
            except OSError:
                pass

    def download_aspera(self, user, host, silent=False):
        """Download file with Aspera Connect.

        For details see the documentation ov Aspera Connect

        Args:
            user (:obj:`str`): FTP user.
            host (:obj:`str`): FTP host. Defaults to "ftp-trace.ncbi.nlm.nih.gov".
        """
        aspera_home = os.environ.get("ASPERA_HOME", None)
        if not aspera_home:
            raise ValueError("environment variable $ASPERA_HOME not set")
        if not os.path.exists(aspera_home):
            raise ValueError(
                "$ASPERA_HOME directory {} does not exist".format(aspera_home)
            )
        ascp = os.path.join(aspera_home, "connect/bin/ascp")
        key = os.path.join(aspera_home, "connect/etc/asperaweb_id_dsa.openssh")
        if not os.path.exists(ascp):
            raise ValueError("could not find ascp binary")
        if not os.path.exists(key):
            raise ValueError("could not find openssh key")

        parsed_url = urlparse(self.url)

        cmd = "{} -i {} -k1 -T -l400m {}@{}:{} {}".format(
            ascp, key, user, host, parsed_url.path, self._temp_file_name
        )
        logger.debug(cmd)
        try:
            pr = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = pr.communicate()
            if not silent:
                logger.debug("Aspera stdout: " + str(stdout))
                logger.debug("Aspera stderr: " + str(stderr))
            if pr.returncode == 0:
                logger.debug(
                    "Moving %s to %s" % (self._temp_file_name, self.destination)
                )
                shutil.move(self._temp_file_name, self.destination)
                logger.debug("Successfully downloaded %s" % self.url)
            else:
                logger.error("Failed to download %s using Aspera Connect" % self.url)
        finally:
            try:
                os.remove(self._temp_file_name)
            except OSError:
                pass

    def _get_filename(self):
        filename = os.path.basename(urlparse(self.url).path).strip(" \n\t.")
        if len(filename) == 0:
            raise Exception("Cannot parse filename from %s" % self.url)
        return filename

    def _download_ftp(self, silent=False):
        parsed_url = urlparse(self.url)
        try:
            ftp = FTP(parsed_url.netloc)
            ftp.login()
            total_size = ftp.size(parsed_url.path)
            if total_size is None:
                total_size = 0
            wrote = list()  # cannot add in the callback, has to be a list
            with open(self._temp_file_name, "wb") as f:
                if silent:

                    def _write(data):
                        f.write(data)
                        wrote.append(len(data))

                    ftp.retrbinary("RETR %s" % parsed_url.path, _write)
                else:
                    with tqdm(
                        total=total_size,
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                        leave=True,
                    ) as pbar:

                        def _write(data):
                            data_length = len(data)
                            pbar.update(data_length)
                            f.write(data)
                            wrote.append(data_length)

                        ftp.retrbinary("RETR %s" % parsed_url.path, _write)
            ftp.quit()
        except Exception:
            try:
                ftp.quit()
                logger.error(
                    "Error when trying to retreive %s." % self.url, exc_info=True
                )
            except Exception:
                logger.error("Error when quiting FTP server.", exc_info=True)

        if total_size != 0:
            if sum(wrote) != total_size:
                raise ValueError(
                    "Downloaded size do not match the expected size for %s" % (self.url)
                )
            else:
                logger.debug("Size validation passed")

    def _download_http(self, silent=False):
        r = requests.get(self.url, stream=True)
        r.raise_for_status()
        # Total size in bytes.
        total_size = int(r.headers.get("content-length", 0))
        logger.debug("Total size: %s" % total_size)
        md5_header = r.headers.get("Content-MD5")
        logger.debug("md5: %s" % str(md5_header))
        chunk_size = 1024
        wrote = 0
        with open(self._temp_file_name, "wb") as f:
            if silent:
                for chunk in r.iter_content(chunk_size):
                    if chunk:
                        f.write(chunk)
                        wrote += len(chunk)
            else:
                with tqdm(
                    total=total_size,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    leave=True,
                ) as pbar:
                    for chunk in r.iter_content(chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
                            wrote += len(chunk)
        if total_size != 0:
            if wrote != total_size:
                raise ValueError(
                    "Downloaded size do not match the expected size for %s" % (self.url)
                )
            else:
                logger.debug("Size validation passed")
        if md5_header:
            logger.debug("Validating MD5 checksum...")
            if md5_header == Downloader.md5sum(self._temp_file_name):
                logger.debug("MD5 checksum passed")
            else:
                raise ValueError("MD5 checksum do NOT passed")

    @staticmethod
    def md5sum(filename, blocksize=8192):
        """Get the MD5 checksum of a file."""
        with open(filename, "rb") as fh:
            m = hashlib.md5()
            while True:
                data = fh.read(blocksize)
                if not data:
                    break
                m.update(data)
        return m.hexdigest()
