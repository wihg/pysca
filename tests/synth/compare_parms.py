from __future__ import print_function
import os, pyfits
from pylab import *

save_plots = 0
fname_plot_scatter = 'compare_params_scatter.pdf'
fname_plot_errors = 'compare_params_errors.pdf'

dpath = '/home/kolja/work/pysca/tests/data/runs/synth'

in_fname_ts = '6761539_scp_syn_n40.fits'
in_fname_parms = 'out_6761539_scp_syn_n40_f14-36_b3_parms.fits'
fname_plot_scatter = '6761539-syn-n40-scatter.pdf'
fname_plot_errors = '6761539-syn-n40-errors.pdf'

#in_fname_ts = '6761539_scp_syn_n40_noi.fits'
#in_fname_parms = 'out_6761539_scp_syn_n40_noi_f14-36_b3_parms.fits'
#fname_plot_scatter = '6761539-syn-n40-noi-scatter.pdf'
#fname_plot_errors = '6761539-syn-n40-noi-errors.pdf'

# load results
p_org = pyfits.getdata(os.path.join(dpath, in_fname_ts), 1).view(recarray)
parms = pyfits.getdata(os.path.join(dpath, in_fname_parms)).view(recarray)
t0 = pyfits.getheader(os.path.join(dpath, in_fname_parms))['t0']

# select parameter set and sort results in respect to frequencies
p = parms[parms.n == max(parms.n)]
p = p[argsort(p.freq)]
p_org = p_org[argsort(p_org.freq)]

# normalize phase
phi = mod(p.phase + p.freq * t0, 1)
phi_org = mod(p_org.phase + p_org.freq * t0, 1)

# create plots
fignr = 0
#figsize = (15, 4.5)
#subplot_params = dict(
    #left=0.06, bottom=0.11, right=0.96, top=0.91, wspace=0.3, hspace=0.2)
figsize = (12, 3.5)
subplot_params = dict(
    left=0.08, bottom=0.13, right=0.96, top=0.92, wspace=0.3, hspace=0.2)

# Scatter plots
fignr += 1; figure(fignr, figsize=figsize); clf()

subplot(131)
title('Frequency')
plot(p_org.freq, p.freq, 'bo')
xlabel('Real frequency [d$^{-1}$]')
ylabel('Measured frequency [d$^{-1}$]')
l = min(p_org.freq.min(), p.freq.min()), max(p_org.freq.max(), p.freq.max())
loffs = (l[1]-l[0])*0.05; l = [l[0]-loffs, l[1]+loffs]
l[1] = 35
xlim(*l); ylim(*l)

subplot(132)
title('Amplitude')
plot(p_org.amp, p.amp, 'go')
xlabel('Real amplitude [ppm]')
ylabel('Measured amplitude [ppm]')
l = min(p_org.amp.min(), p.amp.min()), max(p_org.amp.max(), p.amp.max())
loffs = (l[1]-l[0])*0.05; l = [l[0]-loffs, l[1]+loffs]
l[0] = 0
xlim(*l); ylim(*l)

subplot(133)
title('Phase')
plot(phi_org, phi, 'ro')
xlabel('Real phase')
ylabel('Measured phase')

subplots_adjust(**subplot_params)
draw()

# Error plots
fignr += 1; figure(fignr, figsize=figsize); clf()
subplot(131)
title('Frequency')
plot(p_org.amp, 100*abs(p_org.freq-p.freq)/p_org.freq, 'bo')
xlabel('Amplitude [ppm]')
ylabel('Relative error [%]')
yl = ylim(); ylim(-0.02*yl[1], yl[1])

subplot(132)
title('Amplitude')
plot(p_org.amp, 100*abs(p_org.amp-p.amp)/p_org.amp, 'go')
xlabel('Amplitude [ppm]')
ylabel('Relative error [%]')
yl = ylim(); ylim(-0.02*yl[1], yl[1])

subplot(133)
title('Phase')
plot(p_org.amp, 100*abs(phi_org-phi), 'ro')
xlabel('Amplitude [ppm]')
ylabel('Absolute error [%]')
yl = ylim(); ylim(-0.02*yl[1], yl[1])

subplots_adjust(**subplot_params)
draw()

if save_plots:
    figure(1).savefig(fname_plot_scatter)
    figure(2).savefig(fname_plot_errors)


show()

