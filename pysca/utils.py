from __future__ import absolute_import, division

import numpy as np
from ._fasper import fasper

__author__ = 'Wiebke Herzberg, Kolja Glogowski'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

class ExportDecorator(object):
    def __init__(self, all_list):
        self._all = all_list
    def __call__(self, obj):
        if obj.__name__ not in self._all:
            self._all.append(obj.__name__)
        return obj

__all__ = []
export = ExportDecorator(__all__)

@export
def compute_periodogram(t, a, ofac, hifreq):
    """
    Computes the Lomb-Scargle periodogram with normalized amplitudes.

    Parameters
    ----------
    t : ndarray (1d)
        Time values of the time series.
    a : ndarray (1d)
        Amplitudes of the time series.
    ofac : scalar
        Oversampling factor.
    hifreq : scalar
        Maximum frequency.

    Returns
    -------
    nu : ndarray
        Frequencies of the periodogram.
    psq : ndarray
        Normalized amplitude values of the periodogram.
    """
    t, a = np.atleast_1d(t, a)
    if t.ndim > 1 or a.ndim > 1:
        raise ValueError('Input arrays must be 1d.')
    nu, p, nout, jmax, var = fasper(t, a, float(ofac), hifreq=float(hifreq))

    # Note: amplitudes from fasper are not scaled correctly. To get
    # amplitudes corresponding to the amplitude values from the fit, do
    # the normalization 2.0*sqrt(var*wk2(i)/n) (from kperiod).
    psq = 2.0 * np.sqrt(var * p / len(a))
    return nu, psq

@export
def find_highest_peak(nu, p):
    """
    Find the frequency of the highest peak in the periodogram, using a
    3-point parabolic interpolation.

    Parameters
    ----------
    nu : ndarray (1d)
        Frequency grid of the periodogram.
    p : ndarray (1d)
        Periodogram values.

    Returns
    -------
    nu_peak : scalar
        Interpolated frequency value of the highest peak.
    """
    nu, p = np.atleast_1d(nu, p)
    if nu.ndim > 1 or p.ndim > 1:
        raise ValueError('Input arrays must be 1d.')

    # Get index of highest peak.
    imax = np.argmax(p)

    # Determine the frequency value by parabolic interpolation, which is
    # only possible if the maximum is not at the beginning or end of the
    # periodogram.
    if imax == 0 or imax == p.size - 1:
        nu_peak = p[imax]
    else:
        # Get values around the maximum.
        frq1 = nu[imax-1]
        frq2 = nu[imax]
        frq3 = nu[imax+1]
        y1 = p[imax-1]
        y2 = p[imax]
        y3 = p[imax+1]

        # Parabolic interpolation formula.
        t1 = (y2-y3) * (frq2-frq1)**2 - (y2-y1) * (frq2-frq3)**2
        t2 = (y2-y3) * (frq2-frq1) - (y2-y1) * (frq2-frq3)
        nu_peak = frq2 - 0.5 * t1/t2
    return nu_peak

@export
def median_noise_level(nu, p, nu0, width):
    """
    Computes noise level around a peak using the median of the periodogram.

    Parameters
    ----------
    nu : ndarray (1d)
        Frequency grid of the periodogram.
    p : ndarray (1d)
        Periodogram values.
    nu0 : scalar
        Central frequency of the peak.
    width : scalar
        Area around the central peak to be used for the noise computation.

    Returns
    -------
    noise : scalar
        Noise level.
    """
    nu, p = np.atleast_1d(nu, p)
    if nu.ndim > 1 or p.ndim > 1:
        raise ValueError('Input arrays must be 1d.')
    nu0, width = float(nu0), float(width)
    idx = (nu >= nu0 - 0.5 * width) & (nu <= nu0 + 0.5 * width)
    return np.median(p[idx])
