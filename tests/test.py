import os, pyfits
from pysca import Pysca
from pylab import *

in_dir = 'data/scdata'
#in_fname = 'LC_6761539.fits'
in_fname = 'LC_6443122.fits'
data = pyfits.getdata(os.path.join(in_dir, in_fname))

#in_dir = 'data/lcdata'
#in_fname = 'lc_6761539_tpf_corrconcat-chisq.dat'
#data = loadtxt(os.path.join(in_dir, in_fname))

data = data[isfinite(data[:,1])]
data = data[:36000]
t, a = data[:,0], data[:,1]

p = Pysca(t, a, 0.15, 70.0, 3.0, ofac=6.0)
