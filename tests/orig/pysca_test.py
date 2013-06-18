#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script extracts frequencies, amplitudes and phases. Additionally it
# calculates the signal to noise ratio (S/N) for every frequency. Values you
# need to set : numin, numax, snlimit, amlimit, frqnumlimit, noibox,
# input filename, output root filename.
# The oversampling factor ofac can be left at 10. even with sc data, when
# running on a server. Change down to 1 or 2 when running on a laptop (memory!)
# Special Option: per default, the periodogram gets calculated up to frequency
# 'numax'. If you want to calculate it up to another value, you can set the
# parameter hifreq in the call of the 'lomb_w.fasper()' function manually to
# your desired value.

import time
from numpy import *
import lomb_w
import prewhitefunc
import pyfits

# Limits for the range in which you want the frequencies extracted and the S/N
# and/or amplitude down to which you continue to extract. Set 'snlimit' or
# 'amlimit' to zero if you do not wish to use them.
# Note: the S/N that is calculated in the loop can, and probably will, differ
# from the final S/N that is calculated after the last extraction and fit,
# because all amplitudes are still optmized (and will change) until the end.
numin = 0.15    #0.15 #22.5
numax = 70.0   #22.5 #70.0
snlimit = 5.0
amlimit = 0.0
frqnumlimit = 1

# Box size for noise calculation in c/d. Should be smaller than the range in
# which you look for frequencies (< numax-numin)
noibox = 3.
assert noibox < numax - numin, 'Noise box too big'

#output root name
out_fileroot = '../data/out/pysca_test-'

#load time series
#a = loadtxt('/home/wiebke/data/TPF/6761539/lc/corrconcat/subsets/LC_6761539_tpf_osplit17_part17.dat')
a = pyfits.getdata('../data/scdata/LC_6443122.fits')

a = a[isfinite(a[:,1])] # Check for NaNs that are present in the fits files
n = a.shape[0]
# create a copy to do the prewhitening on
b = a.copy()

#Initialize empty list for frequencies
frqs = []

#Initialize empty list for noise
noise = []

# Set startvalues for fit to None for first run through loop
optpar = None

tstart = time.time()

i = 0 #Counter for periodogram savings
#Loop for a little while (until S/N drops below 'snlimit')
while True:
    tloop = time.time()

    #for ofac an even number must be chosen, otherwise error occurs
    px, py, nout, jmax, var = lomb_w.fasper(b[:,0], b[:,1], 6., hifreq=200)

    # Note: amplitudes from fasper are not scaled correctly. To get amplitudes
    # corresponding to the amplitude values from the fit, do the normalization
    # 2.0*sqrt(var*wk2(i)/n) (from kperiod).
    pysq = 2. * sqrt(var * py / n)
    resp = c_[px, pysq]

    #Write periodogram to fits file...
    if i == 0:
        pyfits.writeto(out_fileroot + 'perd.fits', resp, clobber=True)
    elif i <= 5 or mod(i, 10) == 0:
        pyfits.append(out_fileroot + 'perd.fits', resp, verify=False)

    # Limit the periodogram to selected frequency range
    idx = (px>=numin) & (px<=numax)
    pysqcut = pysq[idx]
    pxcut = px[idx]

    # Calculate noise for last extracted frequency
    if frqs:
        pysq_noibox = pysqcut[(pxcut >= frqs[-1] - noibox/2.) & (pxcut<= frqs[-1] + noibox/2.)]
        # Use median of periodogram values to calculate a noise level
        noi = median(pysq_noibox)
        noise.append(noi)
        # Calculate S/N ratio
        sn = optpar[:,0][-1]/noi
        print 'S/N', sn

        if sn < snlimit or optpar[:,0][-1] < amlimit or i == frqnumlimit:
            break

    # Get new frequency of highest peak in selected range
    frq = pxcut[pysqcut == max(pysqcut)]
    assert frq.size == 1, 'Ambiguity in frequencies detected!'
    print 'Frequency:', frq
    
    # Refine the frequency value by parabolic interpolation.
    # Get index of maximum
    imax = argmax(pysqcut)

    # Handle the case that the maximum is at beginning or end of periodogram
    if imax == 0 or imax == pysqcut.size - 1:
        peak = frq[0]
    else:
        # Get values around the maximum
        frq1 = pxcut[imax - 1]
        frq2 = pxcut[imax]
        frq3 = pxcut[imax + 1]
        y1   = pysqcut[imax - 1]
        y2   = pysqcut[imax]
        y3   = pysqcut[imax + 1]

        # Parabolic interpolation formula
        t1 = (y2-y3) * (frq2-frq1)**2 - (y2-y1) * (frq2-frq3)**2
        t2 = (y2-y3) * (frq2-frq1) - (y2-y1) * (frq2-frq3)
        peak = frq2 - 0.5 * t1/t2

    print '     Peak:', peak

    # Add refined peak frequency to list
    frqs.append(peak)

    # Note: do the prewhitening always with the original time series 'a'!
    b, optpar = prewhitefunc.prewhiten(a, frqs, optpar)
    #print optpar, '\n'
    print len(frqs), 'frequencies fitted.'

    i += 1
    
    now = time.time()
    print 'Elapsed time: %.1fs (loop) / %.1fs (total)\n' % ((now - tloop), (now - tstart))


frqs = array(frqs)
noise = array(noise)
signoi = optpar[:,0]/noise
framph = c_[frqs, optpar, signoi, noise]

# Write extracted frequency data to file
pyfits.writeto(out_fileroot + 'frqs.fits', framph, clobber=True)
# Write residuals of the time series to file
pyfits.writeto(out_fileroot + 'resi.fits', b, clobber=True)

#better use fits files for saving data!
#savetxt('loscaper-newscale.dat', resp, fmt=('%17e', '%17e'))
