#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import re
from glob import glob
from os.path import basename, dirname, join, splitext

from setuptools import find_packages

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read().replace(".. :changelog:", "")

# TODO: Automatically add requirements from requirements.txt
requirements = [
    "numpy>=1.7",
    "pandas>=0.17",
    "requests>=2.21.0",
    "tqdm>=4.31.1",
]

test_requirements = [
    "tox",
]

setup(
    name="GEOparse",
    version="2.0.3",
    description="Python library to access Gene Expression Omnibus Database (GEO)",
    long_description=readme + "\n\n" + history,
    author="Rafal Gumienny",
    author_email="guma44@gmail.com",
    url="https://github.com/guma44/GEOparse",
    packages=find_packages("src"),
    package_dir={"": "src"},
    python_requires=">3.5.0",
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    include_package_data=True,
    scripts=["scripts/geo2fastq"],
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=[
        "GEOparse",
        "GEO",
        "Gene Expression Omnibus",
        "Bioinformatics",
        "Microarray",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    test_suite="tests",
    tests_require=test_requirements,
    project_urls={
        "Documentation": "https://geoparse.readthedocs.io/",
        "Changelog": "https://geoparse.readthedocs.io/en/latest/history.html",
        "Issue Tracker": "https://github.com/guma44/GEOparse/issues",
    },
)
