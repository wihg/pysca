from __future__ import absolute_import, print_function

import os, pyfits, time
from pysca import Pysca
from pysca.io import *
import numpy as np
from numpy import *

##in_fname = 'data/scdata/LC_6761539.fits'
#in_fname = 'data/scdata/LC_6443122.fits'
#in_fname = 'data/lcdata/lc_6761539_tpf_chisq-ppm.dat'
#in_fname = 'test_lc_2.fits'
in_fname = 'test_lc_2_s200.fits'

#n = None
n = 5000
#n = 2000

t, a = read_timeseries(in_fname)
t, a = t[:n], a[:n]

try:
    p = Pysca(t, a, 0.5, 20.0, 3.0, ofac=6)
    for i in range(11):
        print(i)
        p.step()
    res = p.result.view(recarray)
except KeyboardInterrupt:
    print('muh')

freq0 = pyfits.getdata(in_fname, 1).freq[:len(res)]
amp0 = pyfits.getdata(in_fname, 1).amp[:len(res)]
phi0 = pyfits.getdata(in_fname, 1).phi[:len(res)]

ix = []
for i in range(len(res)):
    ix.append(argmin(abs(freq0 - res.freq[i])))
res = res[ix]

print(res.freq)
print(res.amp)
print(res.phase)

x = c_[freq0, freq0-res.freq, amp0-res.amp, phi0-res.phase]
print(x); print()
y = c_[(freq0-res.freq)/freq0, (amp0-res.amp)/amp0, (phi0-res.phase)/phi0]
print(c_[freq0, abs(y) * 100])
