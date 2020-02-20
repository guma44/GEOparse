"""A module that helps to deal with SRA files."""

import glob
import json
import os
import platform
import re
import subprocess as sp

from pandas import DataFrame, concat
from six import iteritems

from . import utils
from .logger import geoparse_logger as logger


class NoSRARelationException(Exception):
    pass


class FastqDumpException(Exception):
    pass


class SRADownloader(object):
    """Manage download RAW data as SRA files.

    The files will be downloaded to the sample directory created ad hoc
    or the directory specified by the parameter. The sample has to come
    from sequencing eg. mRNA-seq, CLIP etc.

    An important parameter is a filetype. By default an SRA
    is accessed by FTP and such file is downloaded. This does not
    require additional libraries. However in order
    to produce FASTA of FASTQ files one would need to use SRA-Toolkit.
    Thus, it is assumed that this library is already installed or it
    will be installed in the near future. One can immediately specify
    the download type to fasta or fastq.
    """

    ALLOWED_FILETYPES = ("sra", "fastq", "fasta")
    FTP_ADDRESS_TPL = (
        "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads"
        "/ByRun/sra/SRR/{range_subdir}/{file_dir}/"
        "{file_dir}.sra"
    )

    def __init__(self, gsm, email, directory="./", **kwargs):
        """Initialize downloader object.

        Args:
            gsm (:class:`GEOparse.GSM`): A GSM object
            email (:obj:`str`): an email (any) - Required by NCBI for access
            directory (:obj:`str`, optional): The directory to which download
                the data. Defaults to "./".
            **kwargs: Arbitrary keyword arguments, see description

        Following  ``**kwargs`` can be passed:
            * filetype - str
                can be sra, fasta, or fastq - for fasta or fastq SRA-Toolkit
                need to be installed
            * aspera - bool
                use Aspera to download samples, defaults to False
            * keep_sra - bool
                keep SRA files after download. Removes SRA file only if the
                selected file type is different than "sra", defaults to False
            * threads - int
                number of threads to use if parallel-fastq-dump is used,
                defaults to 4
            * fastq_dump_options - dict
                pass options to fastq-dump (if used, the options has to be in
                long form eg. --split-files), defaults to::
                    {
                        'split-files': None,
                        'readids': None,
                        'read-filter': 'pass',
                        'dumpbase': None,
                        'gzip': None
                    }

        Raises:
            :obj:`TypeError`: Type to download unknown
            :obj:`TypeError`: Wrong e-mail
        """
        # Unpack arguments
        self.gsm = gsm
        self.email = email
        # Retrieving output directory name

        if platform.system() == "Windows":
            name_regex = r"[\s\*\?\(\),\.\:\%\|\"\<\>]"
        else:
            name_regex = r"[\s\*\?\(\),\.;]"

        self.directory = os.path.abspath(
            os.path.join(
                directory,
                "%s_%s_%s"
                % (
                    "Supp",
                    self.gsm.get_accession(),
                    re.sub(name_regex, "_", self.gsm.metadata["title"][0]),
                ),
            )
        )

        self.filetype = kwargs.get("filetype", "fasta").lower()
        self.aspera = kwargs.get("aspera", False)
        self.keep_sra = kwargs.get("keep_sra", False)
        self.silent = kwargs.get("silent", False)
        self.force = kwargs.get("force", False)
        self.threads = kwargs.get("threads", 4)

        self.fastq_dump_options = {
            "split-files": None,
            "readids": None,
            "read-filter": "pass",
            "dumpbase": None,
            "gzip": None,
        }
        self.fastq_dump_options.update(kwargs.get("fastq_dump_options", {}))

        if self.filetype not in type(self).ALLOWED_FILETYPES:
            raise TypeError(
                "Unknown type to downlod: %s. Allowed filetypes: %s"
                % (self.filetype, type(self).ALLOWED_FILETYPES)
            )

        if not ("@" in email and email != "" and "." in email):
            raise TypeError("Provided e-mail (%s) is invalid" % self.email)
        self._paths_for_download = None

    @property
    def paths_for_download(self):
        """List of URLs available for downloading."""
        if self._paths_for_download is None:
            queries = list()
            try:
                for sra in self.gsm.relations["SRA"]:
                    query = sra.split("=")[-1]
                    if "SRX" not in query:
                        raise ValueError(
                            "Sample looks like it is not an SRA: %s" % query
                        )
                    logger.info("Query: %s" % query)
                    queries.append(query)
            except KeyError:
                raise NoSRARelationException(
                    "No relation called SRA for %s" % self.gsm.get_accession()
                )

            # Construction of DataFrame df with paths to download
            df = DataFrame(columns=["download_path"])
            for query in queries:
                # retrieve IDs for given SRX
                answer = utils.esearch(db="sra", term=query, email=self.email)
                ids = answer["esearchresult"]["idlist"]
                if len(ids) != 1:
                    raise ValueError("There should be one and only one ID per SRX")

                results = utils.efetch(
                    db="sra", id=ids[0], rettype="runinfo", email=self.email
                )
                try:
                    df_tmp = DataFrame(
                        [i.split(",") for i in results.split("\n") if i != ""][1:],
                        columns=[i.split(",") for i in results.split("\n") if i != ""][
                            0
                        ],
                    )
                except IndexError:
                    logger.error(
                        (
                            "SRA is empty (ID: %s, query: %s). "
                            "Check if it is publicly available."
                        )
                        % (ids[0], query)
                    )
                    continue

                # check it first
                try:
                    df_tmp["download_path"]
                except KeyError as e:
                    logger.error("KeyError: " + str(e) + "\n")
                    logger.error(str(results) + "\n")

                df = concat([df, df_tmp], sort=True)
            self._paths_for_download = [path for path in df["download_path"]]
        return self._paths_for_download

    def download(self):
        """Download SRA files.

        Returns:
            :obj:`list` of :obj:`str`: List of downloaded files.
        """
        self.downloaded_paths = list()
        for path in self.paths_for_download:
            downloaded_path = list()
            utils.mkdir_p(os.path.abspath(self.directory))

            sra_run = path.split("/")[-1]
            logger.info("Analysing %s" % sra_run)
            if self.aspera:
                url = type(self).FTP_ADDRESS_TPL.format(
                    range_subdir=sra_run[:6], file_dir=sra_run
                )
            else:
                url = path
            logger.debug("URL: %s", url)
            filepath = os.path.abspath(os.path.join(self.directory, "%s.sra" % sra_run))
            utils.download_from_url(
                url, filepath, aspera=self.aspera, silent=self.silent, force=self.force
            )
            logger.info("Successfully downloaded %s" % filepath)

            if self.filetype in ("fasta", "fastq"):
                if utils.which("fastq-dump") is None:
                    logger.error("fastq-dump command not found")
                ftype = ""
                if self.filetype == "fasta":
                    ftype = " --fasta "
                cmd = "fastq-dump"
                if utils.which("parallel-fastq-dump") is None:
                    cmd += " %s --outdir %s %s"
                else:
                    logger.debug("Using parallel fastq-dump")
                    cmd = " parallel-fastq-dump --threads %s"
                    cmd = cmd % self.threads
                    cmd += " %s --outdir %s -s %s"
                cmd = cmd % (ftype, self.directory, filepath)

                for fqoption, fqvalue in iteritems(self.fastq_dump_options):
                    if fqvalue:
                        cmd += " --%s %s" % (fqoption, fqvalue)
                    elif fqvalue is None:
                        cmd += " --%s" % fqoption
                logger.debug(cmd)
                process = sp.Popen(
                    cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True, bufsize=9999999
                )
                logger.info(
                    "Converting to %s/%s*.%s.gz\n"
                    % (self.directory, sra_run, self.filetype)
                )
                stdout, stderr = process.communicate()
                if process.returncode != 0:
                    raise FastqDumpException("fastq-dump failed: %s" % stderr)
                downloaded_path = glob.glob(
                    os.path.join(self.directory, "%s*.%s.gz" % (sra_run, self.filetype))
                )

            elif self.filetype == "sra":
                downloaded_path = glob.glob(
                    os.path.join(self.directory, "%s*.%s" % (sra_run, self.filetype))
                )

            else:
                downloaded_path = glob.glob(
                    os.path.join(self.directory, "%s*" % sra_run)
                )
                logger.error("Filetype %s not supported." % self.filetype)

            if not self.keep_sra and self.filetype != "sra":
                # Delete sra file
                os.unlink(filepath)

            self.downloaded_paths += downloaded_path
        return self.downloaded_paths
