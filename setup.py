#!/usr/bin/env python

import os, re

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

NAME = 'pysca'
README = open('README.txt').readlines()
DESCRIPTION = README[0]
LONG_DESCRIPTION = '\n'.join(README[2:])
AUTHOR = 'Wiebke Herzberg, Kolja Glogowski'
AUTHOR_EMAIL = '"Kolja Glogowski" <kolja@kis.uni-freiburg.de>'
URL = 'http://pypi.python.org/pypi/pysca'
LICENSE = 'MIT'

# Read version string from pysca/__init__.py without importing the module to
# prevent an import error in case numpy, scipy or pyfits are not installed yet.
m = re.search(r"""^\s*__version__\s*=\s*["'](.+)["']\s*$""",
              open('pysca/__init__.py').read(), re.MULTILINE)
VERSION = m.group(1)

# For distutils builds: Make sure the MANIFEST file ist generated from
# MANIFEST.in by removing the old file.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      url=URL,
      license=LICENSE,
      packages=['pysca'],
      scripts=['bin/pysca'],
      install_requires=['numpy>=1.4.1', 'scipy>=0.7.2', 'pyfits>=2.3.1'],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2.6',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Astronomy'],
      )
