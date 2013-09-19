from __future__ import print_function
import os, pyfits, pysca
from pysca.io import read_timeseries
import numpy.random as rnd
from numpy import *

comp_perds = 1
show_plots = 1
save_synth_ts = 0

#n = 30
n = 40
numin, numax = 14.0, 36.0
hifreq = 72.0
ofac = 6.0
noise_seed = 0

dpath = '/home/kolja/work/pysca/tests/data/runs/synth'

in_fname_ts = '6761539_scp.fits'
in_fname_params = 'out_6761539_scp_f14-36_b3_parms.fits'

out_fname_ts = '6761539_scp_syn_n%d.fits' % n
out_fname_ts_noise = '6761539_scp_syn_n%d_noi.fits' % n

t, a_org = read_timeseries(os.path.join(dpath, in_fname_ts))
#t, a_org = t[:10000], a_org[:10000]
params = pyfits.getdata(os.path.join(dpath, in_fname_params)).view(recarray)

p = params[params.n == n]
a = pysca.harmfunc(t, p.freq, p.amp, p.phase)

noise = sqrt(var(a_org)-var(a))
rnd.seed(noise_seed)
an = a + rnd.normal(0, noise, len(a))
print('Noise level:', noise)

if comp_perds:
    print('Computing original periodogram...')
    nu, psq_org = pysca.compute_periodogram(t, a_org, ofac, hifreq)
    print('Computing clean synthetic periodogram...')
    nu, psq = pysca.compute_periodogram(t, a, ofac, hifreq)
    print('Computing noisy synthetic periodogram...')
    nu, psqn = pysca.compute_periodogram(t, an, ofac, hifreq)

if save_synth_ts:
    print('Writing synthetic time series...')
    pyfits.writeto(os.path.join(dpath, out_fname_ts), c_[t, a], clobber=True)
    pyfits.append(os.path.join(dpath, out_fname_ts), p)
    hdu_noise = pyfits.PrimaryHDU(c_[t, an])
    hdu_noise.header['noise'] = noise
    hdu_noise.writeto(os.path.join(dpath, out_fname_ts_noise), clobber=True)
    pyfits.append(os.path.join(dpath, out_fname_ts_noise), p)

if show_plots and comp_perds:
    from pylab import *
    xl = numin, numax
    yl = 0, 2500
    figsize = 14, 4

    figure(1, figsize=figsize); clf()
    title('Real data')
    plot(nu, psq_org)
    xlim(*xl); ylim(*yl)
    xlabel('Frequency [d$^{-1}$]')
    ylabel('Amplitude [ppm]')
    draw(); show()

    figure(2, figsize=figsize); clf()
    title('Sythetic data, no noise')
    plot(nu, psq)
    xlim(*xl); ylim(*yl)
    xlabel('Frequency [d$^{-1}$]')
    ylabel('Amplitude [ppm]')
    draw(); show()

    figure(3, figsize=figsize); clf()
    title('Sythetic data, noise level: %g' % noise)
    plot(nu, psqn)
    xlim(*xl); ylim(*yl)
    xlabel('Frequency [d$^{-1}$]')
    ylabel('Amplitude [ppm]')
    draw(); show()
