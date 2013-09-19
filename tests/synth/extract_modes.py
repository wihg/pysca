import os, pyfits
from pysca import Pysca
from pysca.io import read_timeseries
from numpy import *

numin, numax = 14.0, 36.0
hifreq = 72.0
noibox = 3.0
ofac = 6.0

dpath = '/home/kolja/work/pysca/tests/data/runs/synth'

#n = 50
#in_fname_ts = '6761539_scp.fits'
#out_fname_base = 'out_6761539_scp_f%g-%g_b%g' % (numin, numax, noibox)

n = 30
in_fname_ts = '6761539_scp_syn_n%d.fits' % n
out_fname_base = 'out_6761539_scp_syn_n%d_f%g-%g_b%g' % (n, numin, numax, noibox)
#in_fname_ts = '6761539_scp_syn_n%d_noi.fits' % n
#out_fname_base = 'out_6761539_scp_syn_n%d_noi_f%g-%g_b%g' % (n, numin, numax, noibox)

out_fname_parms = out_fname_base + '_parms.fits'
out_fname_perds = out_fname_base + '_perds.fits'

t, a = read_timeseries(os.path.join(dpath, in_fname_ts))
#t, a = t[:10000], a[:10000]
p = Pysca(t, a, numin, numax, noibox, ofac=ofac, hifreq=hifreq, verbose=1)

perds = zeros((n+1, len(p.nu)))
perds[0] = p.orig_periodogram
reslist = []
for i in range(1,n+1):
    p.run()
    perds[i] = p.next_periodogram
    for rj in p.result:
        reslist.append((i,) + tuple(rj))
results = rec.fromrecords(reslist, dtype=[('n', 'i4')] + p.result.dtype.descr)

pyfits.writeto(os.path.join(dpath, out_fname_parms), results, clobber=True)
f = pyfits.open(os.path.join(dpath, out_fname_parms), mode='update')
f[0].header['t0'] = (t[0], 'start time of the input time series')
f.close()

pyfits.writeto(os.path.join(dpath, out_fname_perds), p.nu, clobber=True)
pyfits.append(os.path.join(dpath, out_fname_perds), perds)
