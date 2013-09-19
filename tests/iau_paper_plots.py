import os, pyfits
from pylab import *

save_plots = 1

in_dir = '.'
in_fname_mask = 'iau_6443122_sc_f1.6-3.1_b3.dat.perd_%04d.fits'
out_fname = 'iau-6443122-sc.eps'
tstr = 'KIC 6443122, short cadence'

figsize = (10, 4)
#figsize = (10, 3.4)
subplots_layout = dict(
    left=0.075, right=0.97, bottom=0.12, top=0.93, wspace=0.04, hspace=0.10)
#subplots_layout = dict(
    #left=0.075, right=0.97, bottom=0.12, top=0.93, wspace=0.04, hspace=0.10)
xl, yl = (1.57, 3.03), (0, 1450)
tstr = 'KIC 6443122, short cadence'

perds = []
for i in range(4):
    perds.append(pyfits.getdata(os.path.join(in_dir, in_fname_mask % i)))
perds = array(perds)
nu, p = perds[:,:,0], perds[:,:,1]

fig = figure(1, figsize=figsize); clf()
for i in range(4):
    ax = fig.add_subplot(2, 2, i+1)
    if i != 0:
        ax.plot(nu[0], p[0], 'gray')
    ax.plot(nu[i], p[i], 'black')
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    if i == 0 or i == 2:
        ax.set_ylabel('Amplitude [ppm]')
        #pass
    else:
        ax.set_yticklabels([])
    if i == 2 or i == 3:
        ax.set_xlabel('Frequency [d$^{-1}$]')
        #pass
    else:
        ax.set_xticklabels([])

#fig.text(0.5, 0.02, 'Frequency [d$^{-1}$]', ha='center', va='center')
#fig.text(0.015, 0.5, 'Amplitude [ppm]',
         #ha='center', va='center', rotation='vertical')

suptitle(tstr, fontsize='medium')
subplots_adjust(**subplots_layout)
draw()

if save_plots:
    fig.savefig(out_fname)

show()
