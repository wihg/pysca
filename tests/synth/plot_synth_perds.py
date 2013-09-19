from __future__ import print_function
import os, pyfits
from pylab import *

save_plots = 1

n = 40
dpath = '/home/kolja/work/pysca/tests/data/runs/synth'

in_fname_perds = 'out_6761539_scp_f14-36_b3_perds.fits'
in_fname_parms = 'out_6761539_scp_f14-36_b3_parms.fits'
in_fname_perds_syn = 'out_6761539_scp_syn_n40_f14-36_b3_perds.fits'
in_fname_parms_syn = 'out_6761539_scp_syn_n40_f14-36_b3_parms.fits'
in_fname_perds_syn_noi = 'out_6761539_scp_syn_n40_noi_f14-36_b3_perds.fits'
in_fname_parms_syn_noi = 'out_6761539_scp_syn_n40_noi_f14-36_b3_parms.fits'

fname_plot_ps_org = '6761539-syn-ps-org.pdf'
fname_plot_ps = '6761539-syn-ps.pdf'
fname_plot_ps_noi = '6761539-syn-ps-noi.pdf'
fname_plot_ps_log = '6761539-syn-ps-log.pdf'
fname_plot_ps_log_pw = '6761539-syn-ps-log-pw.pdf'

nu_org = pyfits.getdata(os.path.join(dpath, in_fname_perds), 0)
a_org = pyfits.getdata(os.path.join(dpath, in_fname_perds), 1)
parms = pyfits.getdata(os.path.join(dpath, in_fname_parms))
p = parms[parms.n == n]

nu = pyfits.getdata(os.path.join(dpath, in_fname_perds_syn), 0)
a = pyfits.getdata(os.path.join(dpath, in_fname_perds_syn), 1)
nu_noi = pyfits.getdata(os.path.join(dpath, in_fname_perds_syn_noi), 0)
a_noi = pyfits.getdata(os.path.join(dpath, in_fname_perds_syn_noi), 1)
assert all(nu == nu_org) and all(nu == nu_noi)

# create plots
fignr = 0

figsize = (12, 3.5)
subplot_params = dict(left=0.06, bottom=0.15, right=0.96, top=0.96)

figsize2 = (12, 8)
subplot_params2 = dict(
    left=0.07, bottom=0.08, right=0.96, top=0.95, wspace=0.07, hspace=0.08)

fignr += 1; figa = figure(fignr, figsize=figsize); clf()
ax = figa.add_subplot(111)
ax.plot(nu, a_org[0], label='KIC 6761539')
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(-10, 2100)
ax.set_xlabel('Frequency [d$^{-1}$]')
ax.set_ylabel('Amplitude [ppm]')
ax.legend(fontsize='medium')
subplots_adjust(**subplot_params)
draw()

fignr += 1; figb = figure(fignr, figsize=figsize); clf()
ax = figb.add_subplot(111)
ax.plot(nu, a[0], label='synthetic time series, no noise')
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(-10, 2100)
ax.set_xlabel('Frequency [d$^{-1}$]')
ax.set_ylabel('Amplitude [ppm]')
ax.legend(fontsize='medium')
subplots_adjust(**subplot_params)
draw()

fignr += 1; figc = figure(fignr, figsize=figsize); clf()
ax = figc.add_subplot(111)
ax.plot(nu, a_noi[0], label='synthetic time series, gaussian noise')
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(-10, 2100)
ax.set_xlabel('Frequency [d$^{-1}$]')
ax.set_ylabel('Amplitude [ppm]')
ax.legend(fontsize='medium')
subplots_adjust(**subplot_params)
draw()


fignr += 1; figd = figure(fignr, figsize=figsize2); clf()

ax = figd.add_subplot(311)
ax.plot(nu, a_org[0], 'b', label='KIC 6761539')
ax.set_xticklabels([])
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(1e-1, 2e4)
ax.semilogy()
ax.legend(fontsize='medium')

ax = figd.add_subplot(312)
ax.plot(nu, a[0], 'b', label='synthetic time series, no noise')
ax.set_xticklabels([])
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(1e-1, 2e4)
ax.semilogy()
ax.legend(fontsize='medium')
ax.set_ylabel('Amplitude [ppm]')

ax = figd.add_subplot(313)
ax.plot(nu, a_noi[0], 'b', label='synthetic time series, gaussian noise')
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(1e-1, 2e4)
ax.semilogy()
ax.legend(fontsize='medium')
ax.set_xlabel('Frequency [d$^{-1}$]')

subplots_adjust(**subplot_params2)
draw()


fignr += 1; fige = figure(fignr, figsize=figsize2); clf()

ax = fige.add_subplot(311)
ax.plot(nu, a_org[0], 'r', label='KIC 6761539')
ax.plot(nu, a_org[40], 'b')
ax.set_xticklabels([])
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(1e-1, 2e4)
ax.semilogy()
ax.legend(fontsize='medium')

ax = fige.add_subplot(312)
ax.plot(nu, a[0], 'r', label='synthetic time series, no noise')
ax.plot(nu, a[40], 'b')
ax.set_xticklabels([])
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(1e-1, 2e4)
ax.semilogy()
ax.legend(fontsize='medium')
ax.set_ylabel('Amplitude [ppm]')

ax = fige.add_subplot(313)
ax.plot(nu, a_noi[0], 'r', label='synthetic time series, gaussian noise')
ax.plot(nu, a_noi[40], 'b')
xl, yl = ax.set_xlim(14, 36), ax.set_ylim(1e-1, 2e4)
ax.semilogy()
ax.legend(fontsize='medium')
ax.set_xlabel('Frequency [d$^{-1}$]')

subplots_adjust(**subplot_params2)
draw()


if save_plots:
    figa.savefig(fname_plot_ps_org)
    figb.savefig(fname_plot_ps)
    figc.savefig(fname_plot_ps_noi)
    figd.savefig(fname_plot_ps_log)
    fige.savefig(fname_plot_ps_log_pw)

show()
