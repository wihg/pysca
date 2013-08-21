from __future__ import absolute_import, print_function

import os, sys, textwrap
import numpy as np
from . import __version__, Pysca
from .utils import VerbosePrinter
from .io import *

try:
    import argparse
except ImportError:
    from . import _argparse as argparse

__author__ = 'Kolja Glogowski, Wiebke Herzberg'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

DESCRIPTION = 'Pysca is a tool for automated mode parameter extraction from ' \
            + 'photometric time series of heat-driven pulsators.'
VERSION_STRING = 'Pysca ' + __version__

def gt_validator(s, x, t):
    x, value = t(x), t(s)
    if value <= x:
        msg = 'value must be greater than %r' % x
        raise argparse.ArgumentTypeError(msg)
    return value

def ge_validator(s, x, t):
    x, value = t(x), t(s)
    if value < x:
        msg = 'value must be at least %r' % x
        raise argparse.ArgumentTypeError(msg)
    return value

class Application(object):
    def __init__(self, argv=None, name=None):
        self.argv = argv
        self.name = name
        try:
            self._parse_args()
        except SystemExit:
            print()
            raise

    def _parse_args(self):
        parser = argparse.ArgumentParser(
            prog=self.name,
            usage=None,
            #description=DESCRIPTION,
            #epilog=DESCRIPTION,
            add_help=False,
            formatter_class = lambda prog: \
                argparse.HelpFormatter(prog, max_help_position=30))

        g = parser.add_argument_group('Required arguments')
        g.add_argument('tsfile', metavar='TSFILE',
                       help='The file containing the time series. '
                           +'Supported are ASCII and FITS files with two '
                           +'columns, where the first column contains the '
                           +'time and the second column the amplitude.')
        g.add_argument('-f', '--freq', nargs=2, metavar=('MIN','MAX'),
                       type=lambda s: ge_validator(s, 0, float), required=True,
                       help='Frequency range used for peak finding')
        g.add_argument('-b', '--noibox', metavar='NBOX', required=True,
                       type=lambda s: ge_validator(s, 0, float),
                       help='Frequency range used for noise calculation; '
                           +'can be set to 0 if no noise should be calculated')
        g.add_argument('-o', '--outfile', metavar='FILE', required=True,
                       help='Output file for the resulting mode parameters')

        g = parser.add_argument_group('Termination conditions')
        g.add_argument('-n', '--num',
                       type=lambda s: gt_validator(s, 0, int),
                       help='Maximum number of iterations')
        g.add_argument('-N', '--snr',
                       type=lambda s: gt_validator(s, 0, float),
                       help='SNR termination condition')
        g.add_argument('-a', '--amp',
                       type=lambda s: gt_validator(s, 0, float),
                       help='Amplitude termination condition')

        g = parser.add_argument_group('Optional arguments')
        g.add_argument('-s', '--ofac', default=6.0,
                       type=lambda s: ge_validator(s, 2.5, float),
                       help='Oversampling factor used for the Lomb-Scargle '
                           +'periodogram (default: %(default)s)')
        g.add_argument('-i', '--infile', metavar='FILE',
                       help='Read mode parameters file')
        g.add_argument('-p', '--perd', action='store_true',
                       help='Write original and last periodogram to file')
        g.add_argument('-t', '--outfmt', metavar='FMT', default='ascii',
                       choices=['ascii', 'fits', 'plainfits'],
                       help='Output format of the parameter file; '
                           +'supported formats: %(choices)s '
                           +'(default: %(default)s)')
        g.add_argument('-v', '--verbose', action='count', default=1,
                       help='Verbose text output')
        g.add_argument('-h', '--help', action='help',
                       help='Show this help message and exit')
        g.add_argument('--version', action='version', version=VERSION_STRING,
                       help='Show program version and exit')

        args = parser.parse_args(self.argv)
        if args.num == None and args.snr == None and args.amp == None:
            msg = '\n'.join(textwrap.wrap('at least one termination '
                +'condition must be specified using arguments -n/--num, '
                +'-N/--snr and/or -a/--amp'))
            parser.error(msg)
        if args.verbose >= 4:
            print(args)
        self.args = args

    def run(self):
        args = self.args
        vprint = VerbosePrinter(args.verbose)

        # Read time series file.
        vprint("Loading time series '%s'..." % args.tsfile, v=2)
        try:
            t, a = read_timeseries(args.tsfile)
        except (IOError, ValueError) as ex:
            print('Error reading time series: %s.\n' % ex, file=sys.stderr)
            return 1

        # Read parameter file, if specified.
        if args.infile != None:
            vprint("Loading mode parameters '%s'..." % args.infile, v=2)
            try:
                inpar = read_params(
                    args.infile, add_noise_cols=bool(args.noibox))

                # Check if the file looks like a time series.
                if len(inpar) > 1000 and np.all(np.isnan(inpar['phase'])):
                    print('Error: The specified parameter file looks like a '
                          +'time series. Aborting.\n', file=sys.stderr)
                    return 1
            except (IOError, ValueError) as ex:
                print('Error reading parameter file: %s.\n' % ex,
                      file=sys.stderr)
                return 1

        # Initialize Pysca
        vprint('Extracting mode parameters...', v=2)
        p = Pysca(t, a, args.freq[0], args.freq[1], args.noibox,
                  ofac=args.ofac, verbose=args.verbose)
        p.run(args.num)
        print(p.result)

        return 0

if __name__ == '__main__':
    app = Application()
    sys.exit(app.run())
