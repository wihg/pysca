from __future__ import absolute_import, division

import numpy as np
from scipy.optimize import leastsq
from .utils import ExportDecorator

__author__ = 'Wiebke Herzberg, Kolja Glogowski'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

__all__ = []
export = ExportDecorator(__all__)

@export
def harmfunc(t, nu, amp, phi):
    """
    Harmonic model function.

    Parameters
    ----------
    t : ndarray (1d)
        Time steps.
    nu : ndarray (1d)
        Frequencies.
    amp : ndarray (1d)
        Amplitudes.
    phi : ndarray (1d)
        Phases.

    Returns
    -------
    res : ndarray
        Time series of harmonic funtions with frequencies nu, amplitudes amp
        and phases phi.
    """
    t, nu, amp, phi = np.atleast_1d(t, nu, amp, phi)
    if t.ndim > 1 or nu.ndim > 1 or amp.ndim > 1 or phi.ndim > 1:
        raise ValueError('Input arrays must be 1d')
    n = len(nu)
    res = np.zeros(len(t))
    for i in range(n):
        res += amp[i] * np.sin(2 * np.pi * (nu[i] * t + phi[i]))
    return res

def _minfunc(amph, xdat, ydat, nu):
    """
    Residue function to be minimized.

    Parameters
    ----------
    amph : ndarray
        Array containing the amplitudes and phases.
    xdat : ndarray
        Data grid.
    ydat : ndarray
        Data values.
    nu : ndarray
        Frequencies.

    Returns
    -------
    res : ndarray
        Residues between data and model.
    """
    n = len(amph) // 2
    return ydat - harmfunc(xdat, nu, amph[:n], amph[n:])

def _dfunc(amph, xdat, ydat, nu):
    """
    Derivative of minfunc (Jacobian matrix).

    Parameters
    ----------
    amph : ndarray
        Array containing the amplitudes and phases.
    xdat : ndarray
        Data grid.
    ydat : ndarray
        Data values.
    nu : ndarray
        Frequencies.

    Returns
    -------
        The Jacobian matrix of minfunc.
    """
    n = len(amph) // 2
    am, ph = amph[:n], amph[n:]
    res = np.zeros((len(amph), len(xdat)))
    for i in range(n):
        res[i] = -np.sin(2*np.pi * (nu[i]*xdat + ph[i]))
    for i in range(n):
        res[n+i] = -am[i] * np.cos(2*np.pi * (nu[i]*xdat + ph[i])) * 2*np.pi
    return res

def fit_timeseries(t, a, freq, amp0=None, phi0=None, **kwargs):
    """
    Fit a time series using a harmonic function with a fixed set of
    frequencies to determine corresponding amplitudes and phases.

    Parameters
    ----------
    t : array_like
    a : array_like
    freq : array_like
    amp0 : None or array_like
    phi0 : None or array_like
    kwargs : dict
        Additional keyword arguments which are passed directly to the
        scipy.leastsq() function.

    Returns
    -------
    amp : ndarray
        Estimated amplitudes from the fit.
    phi : ndarray
        Estimated (normalized) phases from the fit.
    ok : bool
        True if a solution was found; False otherwise.
    misc_output : tuple (cov_x, infodict, mesg, ier)
        Tuple containing additional results returned by the fit routine. See
        scipy.leastsq() for more informations.
    """
    t, a, freq = np.atleast_1d(t, a, freq)
    amp0 = np.ones_like(freq) if amp0 == None else np.atleast_1d(amp0)
    phi0 = 0.5 * np.ones_like(freq) if phi0 == None else np.atleast_1d(phi0)
    for ary in [t, a, freq, amp0, phi0]:
        if ary.ndim > 1:
            raise ValueError('Input arrays must be 1d')

    # Fill missing initial values for amplitudes and phases.
    if amp0.size < freq.size:
        amp0 = np.concatenate((amp0, np.ones(freq.size - amp0.size)))
    if phi0.size < freq.size:
        phi0 = np.concatenate((phi0, 0.5 * np.ones(freq.size - phi0.size)))

    # Perform leastsq fit.
    x, cov_x, infodict, mesg, ier = leastsq(
        _minfunc, np.concatenate((amp0, phi0)), args=(t, a, freq),
        Dfun=_dfunc, full_output=1, col_deriv=1, **kwargs)

    # The fit was successful, if ier in [1, 2, 3, 4]
    ok = (ier >= 1) and (ier <= 4)

    # Extract amplitudes and phases from the fit result.
    n = len(x) // 2
    amp, phi = x[:n], x[n:]

    # Normalizing the results: Flip sign for negative amplitudes by adjusting
    # the corresponding phases and limit phases to the intervall [0,1). Using:
    # -a*sin(x) == a*sin(x+pi)  and  sin(x) == sin(x+2*pi)
    idx = (amp < 0)
    amp[idx] *= -1.0
    phi[idx] += 0.5
    phi = np.mod(phi, 1)

    return amp, phi, ok, (cov_x, infodict, mesg, ier)

@export
def prewhiten(t, a, freq, amp, phi):
    """
    Prewhitens a time series using a harmonic function with the given
    frequencies, amplitudes and phases.

    Parameters
    ----------
    t : array_like
    a : array_like
    freq : array_like
    amp : array_like
    phi : array_like

    Returns
    -------
    tp : ndarray
        Prewhitened time series.
    """
    t, a, freq, amp, phi = np.atleast_1d(t, a, freq, amp, phi)
    for ary in [t, a, freq, amp, phi]:
        if ary.ndim > 1:
            raise ValueError('Input arrays must be 1d')
    if len(t) != len(a):
        raise ValueError("Arrays 't' and 'a' must have equal size")
    if len(freq) != len(amp) or len(amp) != len(phi):
        raise ValueError("Arrays 'freq', 'amp' and 'phi' must have equal size")
    return a - harmfunc(t, freq, amp, phi)

def prewhiten_compat(t, a, freqs, startpars=None):
    if startpars == None:
        amp0 = phi0 = None
    else:
        amp0, phi0 = startpars[:,0], startpars[:,1]
    amp, phi, ok, misc = fit_timeseries(t, a, freqs, amp0, phi0)
    return prewhiten(t, a, freqs, amp, phi), np.c_[amp, phi]
