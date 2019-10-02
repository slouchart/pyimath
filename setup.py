import os
import re
from setuptools import setup

PACKAGE_NAME = 'pyimath'

# parse version, author and contact from package/module without importing or evaluating the code
version = ''
author, author_email = '', ''
with open(os.path.join(PACKAGE_NAME, '__init__.py')) as fh:
    for line in fh:
        if all([version, author, author_email]):
            break
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group('version')

        m = re.search(r"^__author__ = '(?P<author>[^:]+):(?P<author_email>[^']+)'$", line)
        if m:
            author, author_email = m.group('author'), m.group('author_email')


setup(
    name=PACKAGE_NAME,
    version=version,
    license='MIT',
    description='Pure Python library for finite field arithmetic and polynomial manipulation',
    long_description_content_type='text/markdown',
    long_description=open(os.path.join(os.path.dirname(__file__),
                                         'README.md')).read(),
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
    packages=[PACKAGE_NAME, 'tests'],
    platforms=['any', 'all'],
    test_suite='tests.test_all',
    python_requires='>=3',
)
