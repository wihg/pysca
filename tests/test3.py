from __future__ import absolute_import, print_function

import os, pyfits, time
from pysca import Pysca
from pysca.io import *
import numpy as np
from numpy import *

##in_fname = 'data/scdata/LC_6761539.fits'
#in_fname = 'data/scdata/LC_6443122.fits'
#in_fname = 'data/lcdata/lc_6761539_tpf_chisq-ppm.dat'
in_fname = 'test_lc_2.fits'
#in_fname = 'test_lc_2_s200.fits'

#n = None
n = 5000
verbose = 1

t, a = read_timeseries(in_fname)
t, a = t[:n], a[:n]

p1 = Pysca(t, a, 0.5, 20.0, 3.0, ofac=6, verbose=verbose)
p1.run(3)
r1 = p1.result

x = r1.copy()
x.amp = nan
x.phase = nan
p2 = Pysca(t, a, 0.5, 20.0, 3.0, ofac=6, params=r1, verbose=verbose)
p2.step()
r2 = p2.result
