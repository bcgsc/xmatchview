#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""


from setuptools import setup


with open('README.md') as readme_file:
    readme = readme_file.read()


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

    packages=['xmatchview'],
    package_dir={'xmatchview': '.'},
    package_data={'xmatchview': ['./test/*']},
    scripts=[
        'xmatchview.py',
        'xmatchview-conifer.py'
    ],
    install_requires=[
        'pil>=1.1.7',
    ],

    license="GNU General Public License v3",
    long_description=readme,
    include_package_data=True,
    keywords='Genome visualization DNA sequence alignments Smith Waterman Alignments cross_match xmatchview',
    name='xmatchview',

    # setup_requires=[],
    # test_suite='test',
    # tests_require=[],

    url='https://github.com/bcgsc/xmatchview',
    version='1.1.1',
    zip_safe=False,
)
