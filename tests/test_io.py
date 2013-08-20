import os, pyfits, time
from pysca import Pysca
from pysca.io import *
import numpy as np
from numpy import *

##in_fname = 'data/scdata/LC_6761539.fits'
in_fname_1 = 'data/scdata/LC_6443122.fits'
in_fname_2 = 'data/lcdata/lc_6761539_tpf_chisq-ppm.dat'
#in_fname = 'test_lc_2.fits'

#data = data[isfinite(data[:,1])]
#data = data[:360000]
#data = data[:38664]
#data = data[:30000]

t1, a1 = read_timeseries(in_fname_1)
t2, a2 = read_timeseries(in_fname_2)

pa = loadtxt('test_params_1.dat')
pax = zeros(pa.shape[0], dtype=[('freq','f8'), ('amp','f8'), ('phase','f8')])
for i, cni in enumerate(pax.dtype.names):
    pax[cni] = pa[:,i]

pb = c_[pa, arange(pa.shape[0]), 0.1*arange(pa.shape[0])]
pbx = zeros(pa.shape[0], dtype=[('freq','f8'), ('amp','f8'), ('phase','f8'),
                                ('noise','f8'), ('snr','f8')])
for i, cni in enumerate(pbx.dtype.names):
    pbx[cni] = pb[:,i]

c1 = pyfits.Card('DATE', time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime()),
                 '[UTC] file creation time')
c2 = pyfits.Card('CREATOR', 'pysca v0.0.1', 'program that created this file')

write_params('out_test_io.dat', pbx, fmt='ascii', add_to_header=[c1,c2],
             clobber=True)
write_params('out_test_io.fits', pa, fmt='fits', add_to_header=[c1,c2],
             clobber=True)
write_params('out_test_io_tbl.fits', pbx, fmt='fits-tbl',
             add_to_header=[c1,c2], clobber=True)
write_params('out_test_io_img.fits', pbx, fmt='fits-img',
             add_to_header=[c1,c2], clobber=True)

#p = Pysca(t, a, 0.5, 20.0, 3.0, ofac=6)
#for i in range(7):
    #print i
    #p.step()

#print p.freqs
#print p.amplitudes
#print p.phases

#freq0 = pyfits.getdata(os.path.join(in_dir, in_fname), 1).freq[:len(p.results)]
#amp0 = pyfits.getdata(os.path.join(in_dir, in_fname), 1).amp[:len(p.results)]
#phi0 = pyfits.getdata(os.path.join(in_dir, in_fname), 1).phi[:len(p.results)]
#x = c_[freq0, freq0-p.freqs, amp0-p.amplitudes, phi0-p.phases]
#print x
