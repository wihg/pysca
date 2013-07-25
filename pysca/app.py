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
            description='Extracts frequencies, amplitudes, phases and SNR from time series.',
            prog=self.name,
            version=pysca.__version__)
        parser.formatter.max_help_position = 30
        parser.add_option('-v', '--verbose', action='count', default=0,
                dest='verbose', help='verbose text output, can be used multiple times')
        self.opts, self.args = parser.parse_args(args=self.argv)
        print self.opts
        print self.args

    def _load_data(self):
        pass

    def _save_data(self):
        pass

    def run(self):
        return 0
