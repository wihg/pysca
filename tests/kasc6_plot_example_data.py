import pyfits
from pylab import *

save_plots = 1

in_fname = 'out_6761539_lc_16.3_16.9_3.fits'
out_fname = '6761539-lc.pdf'
tstr = 'KIC 6761539, long cadence'
xl, yl = (16.285, 16.915), (-5, 600)

#in_fname = 'out_6443122_sc_1.6_3.0_3.fits'
#out_fname = '6443122-sc.pdf'
#tstr = 'KIC 6443122, short cadence'
#xl, yl = (1.58, 3.02), (-14, 1450)

f = pyfits.open(in_fname)
a, h = f[0].data, f[0].header

#figsize = (12.5, 6)
figsize = (11, 5)
fig = figure(1, figsize=figsize); clf()

for i in range(4):
    ax = fig.add_subplot(2, 2, i+1)
    if i != 0:
        ax.plot(a[0], a[1], 'r')
    ax.plot(a[0], a[i+1], 'b')
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    if i == 0 or i == 2:
        ax.set_ylabel('Amplitude [ppm]')
    else:
        ax.set_yticklabels([])
    if i == 2 or i == 3:
        ax.set_xlabel('Frequency [d$^{-1}$]')
    else:
        ax.set_xticklabels([])

suptitle(tstr, fontsize='medium')
#subplots_adjust(left=0.07, right=0.97, bottom=0.08, top=0.95,
#                wspace=0.045, hspace=0.08)
subplots_adjust(left=0.07, right=0.97, bottom=0.095, top=0.94,
                wspace=0.045, hspace=0.095)
draw()

if save_plots:
    fig.savefig(out_fname)

show()
