from __future__ import print_function
import os, sys, pyfits
from numpy import *

random.seed(0)

argv = sys.argv[:]
if len(argv) not in [4, 5]:
    print('usage: ', os.path.basename(argv[0]),
          '<params> <timesteps> <output> [<sigma>]')
    sys.exit()
in_params_fname = argv[1]
in_tsteps_fname = argv[2]
out_fname = argv[3]
sigma = float(argv[4]) if len(argv) == 5 else 0.0

t = pyfits.getdata(in_tsteps_fname)
params = loadtxt(in_params_fname, ndmin=2)
freq, amp, phi = params[:,0], params[:,1], params[:,2]

a = zeros(len(t))
for i in range(len(freq)):
    a += amp[i] * sin(2*pi * (freq[i]*t + phi[i]))
if sigma > 0:
    a += random.normal(0, sigma, len(t))

hdu = pyfits.PrimaryHDU(c_[t, a])
hdu_params = pyfits.BinTableHDU(
    rec.fromarrays([freq, amp, phi, ones(len(freq)) * sigma],
                   names=['freq', 'amp', 'phi', 'sigma']),
    name='PARAMS')
pyfits.HDUList([hdu, hdu_params]).writeto(out_fname)
