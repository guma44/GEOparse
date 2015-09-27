#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
                'numpy>=1.7',
                'pandas>=0.15',
                ]

test_requirements = [
                'tox',
]

setup(
    name='GEOparse',
    version='0.1.4',
    description="Python library to access Gene Expression Omnibus Database (GEO)",
    long_description=readme + '\n\n' + history,
    author="Rafal Gumienny",
    author_email='guma44@gmail.com',
    url='https://github.com/guma44/GEOparse',
    packages=[
        'GEOparse',
    ],
    package_dir={'GEOparse':
                 'GEOparse'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=['GEOparse', 'GEO', 'Gene Expression Omnibus', 'Bioinformatics', 'Microarray'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        # "Programming Language :: Python :: 2",
        # 'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.3',
        # 'Programming Language :: Python :: 3.4',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
