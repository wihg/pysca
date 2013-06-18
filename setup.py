#!/usr/bin/env python

from distutils.core import setup
from pysca import __version__ as VERSION

NAME = 'pysca'
README = open('README.txt').readlines()
DESCRIPTION = README[0]
LONG_DESCRIPTION = '\n'.join(README[2:])
AUTHOR = 'Wiebke Herzberg, Kolja Glogowski'
AUTHOR_EMAIL = '"Kolja Glogowski" <kolja@kis.uni-freiburg.de>'
URL = 'http://pypi.python.org/pypi/pysca'
LICENSE = 'MIT'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      url=URL,
      license=LICENSE,
      packages=['pysca'],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2.6',
          'Topic :: Scientific/Engineering :: Astronomy'],
      )
