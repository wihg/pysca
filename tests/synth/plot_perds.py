from __future__ import print_function
import os, pyfits
from pylab import *

save_plots = 1
fname_plot_full = 'plot_perds_full.pdf'
fname_plot_pmodes = 'plot_perds_pmodes.pdf'
fname_plot_pmodes_hi = 'plot_perds_pmodes_hi.pdf'
fname_plot_peakex = 'plot_perds_peakex.pdf'
fname_plot_pw_lin = 'plot_perds_pw_lin.pdf'
fname_plot_pw_log = 'plot_perds_pw_log.pdf'

n = 40
dpath = '/home/kolja/work/pysca/tests/data/runs/synth'

in_fname_perds = 'out_6761539_scp_f14-36_b3_perds.fits'
in_fname_parms = 'out_6761539_scp_f14-36_b3_parms.fits'
fname_plot_full = '6761539-ps-full.pdf'
fname_plot_pmodes = '6761539-ps-pmodes.pdf'
fname_plot_pmodes_hi = '6761539-ps-pmodes-hi.pdf'
fname_plot_peakex = '6761539-ps-peakex.pdf'
fname_plot_pw_lin = '6761539-ps-pw-lin.pdf'
fname_plot_pw_log = '6761539-ps-pw-log.pdf'

#in_fname_perds = 'out_6761539_scp_syn_n40_f14-36_b3_perds.fits'
#in_fname_parms = 'out_6761539_scp_syn_n40_f14-36_b3_parms.fits'
#in_fname_perds = 'out_6761539_scp_syn_n40_noi_f14-36_b3_perds.fits'
#in_fname_parms = 'out_6761539_scp_syn_n40_noi_f14-36_b3_parms.fits'

nu = pyfits.getdata(os.path.join(dpath, in_fname_perds), 0)
a = pyfits.getdata(os.path.join(dpath, in_fname_perds), 1)
parms = pyfits.getdata(os.path.join(dpath, in_fname_parms))
p = parms[parms.n == n]
cd2mhz = 1e3 / (24.0 * 3600.0)

# create plots
fignr = 0

#figsize = (15, 4.5)
#subplot_params = dict(left=0.06, bottom=0.11, right=0.96, top=0.90)
figsize = (12, 3.5)
subplot_params = dict(left=0.06, bottom=0.15, right=0.96, top=0.86)

#figsize2 = (15, 9)
figsize2 = (12, 7)
subplot_params2 = dict(
    left=0.07, bottom=0.08, right=0.96, top=0.95, wspace=0.07, hspace=0.08)

figsize3 = (12, 3.5)
subplot_params3 = dict(left=0.06, bottom=0.15, right=0.96, top=0.96)

fignr += 1; figa = figure(fignr, figsize=figsize); clf()
axa = figa.add_subplot(111)
axa.plot(nu, a[0], label='KIC 6761539')
xl, yl = axa.set_xlim(0, 60), axa.set_ylim(0, 2100)
axa.set_xlabel('Frequency [d$^{-1}$]')
axa.set_ylabel('Amplitude [ppm]')
axa.legend(fontsize='medium')
axa2 = axa.twiny()
axa2.set_xlim(xl[0]*cd2mhz, xl[1]*cd2mhz)
axa2.set_xlabel('Frequency [mHz]')
subplots_adjust(**subplot_params)
draw()

fignr += 1; figb = figure(fignr, figsize=figsize); clf()
axb = figb.add_subplot(111)
axb.plot(nu, a[0], label='KIC 6761539, p-modes')
xl, yl = axb.set_xlim(14, 36), axb.set_ylim(0, 2100)
axb.set_xlabel('Frequency [d$^{-1}$]')
axb.set_ylabel('Amplitude [ppm]')
axb.legend(fontsize='medium')
axb2 = axb.twiny()
axb2.set_xlim(xl[0]*cd2mhz, xl[1]*cd2mhz)
axb2.set_xlabel('Frequency [mHz]')
subplots_adjust(**subplot_params)
draw()

fignr += 1; figbh = figure(fignr, figsize=figsize); clf()
numin, numax = 16, 16.93
axbh = figbh.add_subplot(111)
axbh.plot(nu, a[0], label='KIC 6761539, p-modes')
idx = (nu > numin) & (nu < numax)
axbh.plot(nu[idx], a[0][idx], 'r')
xl, yl = axbh.set_xlim(14, 36), axbh.set_ylim(0, 2100)
axbh.set_xlabel('Frequency [d$^{-1}$]')
axbh.set_ylabel('Amplitude [ppm]')
axbh.legend(fontsize='medium')
axb2 = axbh.twiny()
axb2.set_xlim(xl[0]*cd2mhz, xl[1]*cd2mhz)
axb2.set_xlabel('Frequency [mHz]')
subplots_adjust(**subplot_params)
draw()

fignr += 1; figc = figure(fignr, figsize=figsize2); clf()
numin, numax = 16, 16.93
fid = r_[-1, where((p.freq > numin) & (p.freq < numax))[0]]
for i, ni in enumerate(fid+1):
    axci = figc.add_subplot(2, 3, i+1)
    if i != 0:
        vlines(p.freq[fid[1:i+1]], 0, p.amp[fid[1:i+1]], 'r', lw=1.5)
    axci.plot(nu, a[ni], 'b', label='n=%d' % ni)
    #axci.legend(loc='upper left', fontsize='medium')
    axci.set_xlim(numin, numax)
    axci.set_ylim(-1, 750)
    if i not in [0, 3]:
        axci.set_yticklabels([])
    if i in [0, 1, 2]:
        axci.set_xticklabels([])
figc.subplots_adjust(**subplot_params2)
figc.text(0.5, 0.02, 'Frequency [d$^{-1}$]', ha='center', va='center')
figc.text(0.015, 0.5, 'Amplitude [ppm]',
         ha='center', va='center', rotation='vertical')
draw()

fignr += 1; figd = figure(fignr, figsize=figsize3); clf()
axd = figd.add_subplot(111)
axd.plot(nu, a[0], 'r', label='40 extracted peaks')
axd.plot(nu, a[40])
#axd.vlines(p.freq, 0.1, p.amp, 'r')
xl, yl = axd.set_xlim(14, 36), axd.set_ylim(-10, 2100)
axd.set_xlabel('Frequency [d$^{-1}$]')
axd.set_ylabel('Amplitude [ppm]')
axd.legend(fontsize='medium')
subplots_adjust(**subplot_params3)
draw()

fignr += 1; fige = figure(fignr, figsize=figsize3); clf()
axe = fige.add_subplot(111)
axe.plot(nu, a[0], 'r', label='40 extracted peaks')
axe.plot(nu, a[40])
axe.semilogy()
xl, yl = axe.set_xlim(14, 36), axe.set_ylim(1e-1, 2e4)
axe.set_xlabel('Frequency [d$^{-1}$]')
axe.set_ylabel('Amplitude [ppm]')
axe.legend(loc='upper right', fontsize='medium')
subplots_adjust(**subplot_params3)
draw()

if save_plots:
    figa.savefig(fname_plot_full)
    figb.savefig(fname_plot_pmodes)
    figbh.savefig(fname_plot_pmodes_hi)
    figc.savefig(fname_plot_peakex)
    figd.savefig(fname_plot_pw_lin)
    fige.savefig(fname_plot_pw_log)

show()
