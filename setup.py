import os
import re
from setuptools import setup

# parse version from package/module without importing or evaluating the code
version = ''
with open('imath/__init__.py') as fh:
    for line in fh:
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group('version')
            break

setup(
    name='imath',
    version=version,
    license='MIT',
    description='Pure Python Library for finite field arithmetics and polynomial manipulation',
    long_description=open(os.path.join(os.path.dirname(__file__),
                                         'README.md')).read(),
    keywords='python algebra polynomial integers',
    url='https://github.com/slouchart/imath',
    author='Sébastien Louchart',
    author_email='sebastien.louchart@gmail.com',
    classifiers=['License :: OSI Approved :: MIT License',
                   'Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 3.7',
                  ],
    packages=['imath'],
)