#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()


requirements = [
    'pil>=1.1.7',
]

setup_requirements = []

test_requirements = []

setup(
    author="René L. Warren",
    author_email='rwarren@bcgsc.ca',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ],
    description="Cleavage site prediction via de novo assembly",
    scripts=[
        'xmatchview.py',
        'xmatchview-conifer.py'
    ],
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme,
    include_package_data=True,
    keywords='Genome visualization DNA sequence alignments Smith Waterman Alignments cross_match xmatchview',
    name='xmatchview',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/bcgsc/xmatchview',
    version='1.1.1',
    zip_safe=False,
)
