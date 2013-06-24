import os, pyfits
from pysca import Pysca
from pylab import *

show_plots = 1
save_results = 1

#n = 3
#numin, numax, hifreq = 16.3, 16.9, 20.0
#snr_width, ofac = 3.0, 8.0
#in_dir = 'data/lcdata'
#in_fname = 'lc_6761539_tpf_chisq-ppm.dat'
#out_fname = 'out_6761539_lc_%.1f_%.1f_%d.fits' % (numin, numax, n)
#data = loadtxt(os.path.join(in_dir, in_fname))

n = 3
numin, numax, hifreq = 1.6, 3.0, 30.0
snr_width, ofac = 3.0, 8.0
in_dir = 'data/scdata'
in_fname = 'LC_6443122.fits'
out_fname = 'out_6443122_sc_%.1f_%.1f_%d.fits' % (numin, numax, n)
data = pyfits.getdata(os.path.join(in_dir, in_fname))

data = data[isfinite(data[:,1])]
t, a = data[:,0], data[:,1]

per = []
p = Pysca(t, a, numin, numax, snr_width, ofac=ofac, hifreq=hifreq)

print 'computing original lsp...'
per.append(p.orig_periodogram)

for i in range(n):
    print 'extracting freq #%d...' % i
    p.step()
    print p.freqs
    per.append(p.next_periodogram)

if save_results:
    nuper = array([p.nu] + per)
    per_hdu = pyfits.PrimaryHDU(nuper)
    h = per_hdu.header
    h['numin'], h['numax'], h['hifreq'] = numin, numax, hifreq
    h['snr_wid'], h['ofac'] = snr_width, ofac
    res_hdu = pyfits.BinTableHDU(p.results)
    res_hdu.name = 'results'
    hdus = pyfits.HDUList([per_hdu, res_hdu])
    hdus.writeto(out_fname, clobber=True)

if show_plots:
    maxamp = p.amplitudes.max()
    figure(1); clf()
    plot(p.nu, per[0], 'r'); plot(p.nu, per[-1], 'b')
    xlim(numin, numax)
    ylim(-0.02 * maxamp, 1.1 * maxamp)
    draw()
    show()
