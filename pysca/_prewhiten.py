from __future__ import absolute_import

from numpy import *
from scipy.optimize import leastsq

__author__ = 'Wiebke Herzberg'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

# This module can be used to prewhiten a data set with a number of frequencies.
# As input are needed:
# 'data' -> a time series in the form of a two dimensional array, first column
# should be time and second column the flux/magnitude values
# 'freqs' -> a list of frequencies, will be converted into an array internally

# model function
def harmfunc(t, amph, fr):
    # extract amplitudes and phases from amph
    am = amph[:amph.size/2]
    ph = amph[amph.size/2:]

    if fr.size == am.size:
        res = zeros(t.size)
        for j in range(fr.size):
            res += am[j] * sin(2*pi * (fr[j]*t + ph[j]))
        return res
    else:
        print "Amplitudes/phases do not match number of frequencies"
        return None

# function to be minimized
def minfunc(amph, xdat, ydat, fr):
    return ydat - harmfunc(xdat, amph, fr)

# Derivative of minfunc (Jacobi matrix)
def dfunc(amph, xdat, ydat, fr):
    am = amph[:amph.size/2]
    ph = amph[amph.size/2:]
    ds = zeros((amph.size, xdat.size))
    for i in range(amph.size/2):
        ds[i] = -sin(2*pi * (fr[i]*xdat + ph[i]))
    for o in range(amph.size/2):
        ds[amph.size/2 + o] = -am[o] * cos(2*pi * (fr[o]*xdat + ph[o])) * 2*pi
    return ds

def prewhiten(t, a, freqs, startpars=None):
    # Note: if 'startpars' are given, it needs to be a 2-dimensional array,
    # where the first column contains values for amplitudes and the second
    # column contains values for phases. If the array contains not enough values
    # for all frequencies, the given values will be used for the first few
    # frequencies and the rest of the array will be filled with 1. as start
    # values.

    freqs = array(freqs)

    if startpars == None:
        optpar, cov, infodict, mesg, ier = leastsq(minfunc, r_[ones(freqs.size), 0.5*ones(freqs.size)], args=(t, a, freqs), full_output=1, Dfun= dfunc, col_deriv=1)
    else:
        # Make sure start parameter array has the right shape
        assert startpars.shape[1] == 2, 'Start parameter array does not have two columns'
        assert startpars.shape[0] <= freqs.size, 'Start parameter array too big'

        # Pad the array with 1. if it is not large enough (contains only start
        # values for the first few frequencies)
        if startpars.shape[0] < freqs.size:
            pad = ones((freqs.size - startpars.shape[0], 2))
            startpars = r_[startpars, pad]
            #print 'start parameters\n', startpars

        optpar, cov, infodict, mesg, ier = leastsq(minfunc, startpars.T.flatten(), args=(t, a, freqs), full_output=1, Dfun= dfunc, col_deriv=1)

    # Modify amplitudes and phases to get only positive amplitudes as well as phases
    # in the range from 0 to 1.
    # Create copy of amplitudes and phases to keep original values from fit
    optam = copy(optpar[:optpar.size/2])
    optph = copy(optpar[optpar.size/2:])

    # When turning negative amplitudes into positive amplitudes, their corresponding
    # phase has to be shifted about 0.5 (or pi).
    optph[where(optam < 0)] += 0.5

    # Limit the phase to the intervall [0,1)
    normphase = mod(optph, 1)

    # Put positive amplitudes and adjusted phases back together
    normoptpar = r_[abs(optam), normphase]
    # Print number of iterations of the fit
    #print 'Number of iterations:', infodict['nfev']

    # Calculate residuals
    fitcurve = harmfunc(t, normoptpar, freqs)
    residu = a - fitcurve

    return residu, c_[abs(optam), normphase]
