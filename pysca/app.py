from __future__ import absolute_import

import os, sys, pyfits, pysca
from optparse import OptionParser

__author__ = 'Kolja Glogowski'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

class Application(object):
    def __init__(self, argv=None, name=None):
        if argv == None:
            argv = sys.argv
        self.argv = list(argv)
        self.name = name
        self._parse_opts()

    def _parse_opts(self):
        parser = OptionParser(
            usage='usage: %prog [options]',
            description='Blah',
            prog=self.name,
            version=pysca.__version__)
        self.opts, self.args = parser.parse_args(args=self.argv)

    def run(self):
        return 0
