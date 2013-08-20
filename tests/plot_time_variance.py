import os, pyfits
from pysca import Pysca
from pylab import *

#in_dir = 'data/scdata'
##in_fname = 'LC_6761539.fits'
#in_fname = 'LC_6443122.fits'
#data = pyfits.getdata(os.path.join(in_dir, in_fname))

in_dir = 'data/lcdata'
in_fname = 'lc_6761539_tpf_chisq-ppm.dat'
data = loadtxt(os.path.join(in_dir, in_fname))

data = data[isfinite(data[:,1])]
t, a = data[:,0], data[:,1]

dt = t[1:] - t[:-1]
idx = (dt < 1.2 * mean(dt))

figure(1); clf()
clf(); plot((t-t[0])[1:][idx], dt[idx])
draw()
show()
