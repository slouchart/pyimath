# -*- coding: utf-8 -*-

import os
import re
from setuptools import setup, find_packages

PACKAGE_NAME = 'pyimath'
__version__ = '0.1.0'
__author__ = 'Sébastien Louchart:sebastien.louchart@gmail.com'


# parse version, author and contact from package/module without importing or evaluating the code
author, author_email = '', ''

m = re.search(r"^'(?P<author>[^:]+):(?P<author_email>[^']+)'$", __author__)
if m:
    author, author_email = m.group('author'), m.group('author_email')


setup(
    name=PACKAGE_NAME,
    version=__version__,
    license='MIT',
    description='Pure Python library for finite field arithmetic and polynomial manipulation',
    long_description_content_type='text/markdown',
    long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
    keywords='python math arithmetic algebra polynomials number_theory integers',
    url='https://github.com/slouchart/pyimath',
    author='Sébastien Louchart',
    author_email='sebastien.louchart@gmail.com',
    classifiers=['License :: OSI Approved :: MIT License',
                   'Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 3 :: Only',
                   'Intended Audience :: Education',
                   'Operating System :: OS Independent',
                   'Topic :: Scientific/Engineering :: Mathematics',
                  ],
    packages=find_packages(PACKAGE_NAME, exclude=['tests']),
    platforms=['any', 'all'],
    python_requires='>=3',
)
