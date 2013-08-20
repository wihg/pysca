import os, pyfits
from pysca import Pysca
from pylab import *

#in_dir = 'data/scdata'
##in_fname = 'LC_6761539.fits'
#in_fname = 'LC_6443122.fits'
#data = pyfits.getdata(os.path.join(in_dir, in_fname))

#in_dir = 'data/lcdata'
#in_fname = 'lc_6761539_tpf_chisq-ppm.dat'
#data = loadtxt(os.path.join(in_dir, in_fname))

in_dir = '.'
in_fname = 'test_lc_2.fits'
data = pyfits.getdata(os.path.join(in_dir, in_fname))

data = data[isfinite(data[:,1])]
#data = data[:360000]
#data = data[:38664]
#data = data[:30000]
t, a = data[:,0], data[:,1]

#p = Pysca(t, a, 0.15, 70.0, 3.0, ofac=6.0)
#p = Pysca(t, a, 1.5, 3.0, 1.0, ofac=6.0)

## _fasper.__spread__ seems to have problems, if the length of the ts is lower
## than the created array for which the length depend on ofac and hifreq/hifac.
## It seems that an ofac >= 2.1 works for all cases. This needs to be checked!
## -> apparently some time series works with lower ofac values...
# p = Pysca(t, a, 0.5, 20.0, 3.0, ofac=1.0)

p = Pysca(t, a, 0.5, 20.0, 3.0, ofac=10)
for i in range(7):
    print i
    p.step()

print p.freqs
print p.amplitudes
print p.phases

freq0 = pyfits.getdata(os.path.join(in_dir, in_fname), 1).freq[:len(p.results)]
amp0 = pyfits.getdata(os.path.join(in_dir, in_fname), 1).amp[:len(p.results)]
phi0 = pyfits.getdata(os.path.join(in_dir, in_fname), 1).phi[:len(p.results)]
x = c_[freq0, freq0-p.freqs, amp0-p.amplitudes, phi0-p.phases]
print x

#print 100 * abs((p.amplitudes - amp0) / amp0)
#print 100 * abs((p.phases - phi0) / phi0)

#for ri in p.results:
    #print '%.2f  %.2f  %.2f' % (ri.freq, ri.amp, ri.phase)
