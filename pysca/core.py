from __future__ import absolute_import, print_function

import numpy as np
from . import utils, fit
from .utils import ExportDecorator, VerbosePrinter

__author__ = 'Kolja Glogowski, Wiebke Herzberg'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

__all__ = []
export = ExportDecorator(__all__)

@export
class PyscaError(Exception):
    pass

@export
class Pysca(object):
    def __init__(self, t, a, numin, numax, snr_width, ofac=6.0, hifreq=None,
                 verbose=0):
        """
        Pysca(t, a, numin, numax, snr_width, ofac=6.0, hifreq=None)

        Parameters
        ----------
        t : array
            time values of the time series
        a : array
            amplitues of the time series
        numin : float
            minimum frequency value
        numax : float
            maximum frequency value
        snr_width : float
            frequency width around the extracted peak which is used for the
            signal-to-noise calculation
        ofac : float
            oversampling factor used for the Lomb-Scargle periodogram
        hifreq : None
            maximum frequency of the Lomb-Scargle periodogram; default: numax
        verbose : int
            verbose level from 0 (quiet) to 2 (verbose); default 0
        """
        self._t = np.asarray(t, dtype=np.float64)
        self._a = np.asarray(a, dtype=np.float64)
        if self._t.ndim != 1 or self._a.ndim != 1:
            raise ValueError('Input arrays must be 1d')
        elif len(self._t) != len(self._a):
            raise ValueError('Input arrays must have the same sizes')
        self._numin, self._numax = float(numin), float(numax)
        self._snr_width = float(snr_width) if snr_width != None else None
        self._ofac = float(ofac)
        self._hifreq = float(hifreq) if hifreq != None else self._numax
        self._verbose = int(verbose)
        self._vprint = VerbosePrinter(self._verbose)
        self._prev_ts = self._next_ts = None
        self._nu = self._orig_per = self._prev_per = self._next_per = None
        self._freqs = []      # List of extracted frequencies
        self._noise = []      # List of corresponding noise levels
        self._optpar = None   # Start values (am, ph) for the fit (prewhiten)
        self._results = None  # Caching variable for results

    @property
    def numin(self):
        return self._numin

    @property
    def numax(self):
        return self._numax

    @property
    def snr_width(self):
        return self._snr_width

    @property
    def ofac(self):
        return self._ofac

    @property
    def hifreq(self):
        return self._hifreq

    @property
    def t(self):
        return self._t

    @property
    def orig_ts(self):
        return self._a

    @property
    def prev_ts(self):
        return self._prev_ts

    @property
    def next_ts(self):
        return self._next_ts

    @property
    def nu(self):
        if self._nu == None:
            self._nu, self._orig_per = self._calc_periodogram(self._t, self._a)
            self._nuidx = (self._nu >= self.numin) & (self._nu <= self.numax)
        return self._nu

    @property
    def orig_periodogram(self):
        if self._orig_per == None:
            tmp = self.nu
            assert self._nu != None and self._orig_per != None
        return self._orig_per

    @property
    def prev_periodogram(self):
        return self._prev_per

    @property
    def next_periodogram(self):
        return self._next_per

    @property
    def count(self):
        return len(self._freqs)

    @property
    def result(self):
        if self._results == None:
            self._results = np.rec.fromarrays(
                [ self.freqs, self.amplitudes, self.phases, self.noise,
                  self.snr],
                names=['freq', 'amp', 'phase', 'noise', 'snr'])
        return self._results

    @property
    def freqs(self):
        return np.array(self._freqs, dtype=np.float64)

    @property
    def amplitudes(self):
        return self._optpar[:,0] if self._optpar != None else np.array([])

    @property
    def phases(self):
        return self._optpar[:,1] if self._optpar != None else np.array([])

    @property
    def noise(self):
        return np.array(self._noise, dtype=np.float64)

    @property
    def snr(self):
        return self.amplitudes / self.noise

    def _calc_periodogram(self, t, a):
        return utils.compute_periodogram(t, a, self.ofac, self.hifreq)

    def _find_highest_peak(self, nu, per, use_nuidx=True):
        if use_nuidx:
            # Limit the periodogram to selected frequency range.
            nu = nu[self._nuidx]
            per = per[self._nuidx]
        return utils.find_highest_peak(nu, per)

    def step(self):
        if self._next_per == None:
            self._next_per = self.orig_periodogram
            self._next_ts = self.orig_ts

        # Select the next periodogram
        nu, per = self.nu, self.next_periodogram

        # Find frequency of the highest peak in the current periodogram
        freq = self._find_highest_peak(nu, per)
        self._freqs.append(freq)

        # Prewhiten the original time series with all extracted frequencies
        new_ts, new_optpar = fit.prewhiten_compat(
            self.t, self.orig_ts, self._freqs, self._optpar)

        # amplitude of the previously extracted frequency
        amp = new_optpar[-1][0]

        # Compute new periodogram from the prewhitened time series
        new_nu, new_per = self._calc_periodogram(self.t, new_ts)

        # Compute noise for the last extracted frequency using the median
        noise = utils.median_noise_level(nu, new_per, freq, self.snr_width)
        self._noise.append(noise)

        # Update object data
        self._prev_per = self._next_per
        self._prev_ts = self._next_ts
        self._next_per = new_per
        self._next_ts = new_ts
        self._optpar = new_optpar

        # clear results cache
        self._results = None

    def run(self, n=1, amp_limit=None, snr_limit=None):
        """
        Parameters
        ----------
        n : int or None
            number of steps; default: no limit
        amp_limit : float or None
            amplitude limit (termination condition); default: no limit
        snr_limit : float or None
            signal-to-noise limit (termination condition); default: no limit

        Return
        ------
            number of iterations
        """
        n = int(n) if n != None else None
        amp_limit = float(amp_limit) if amp_limit != None else None
        snr_limit = float(snr_limit) if snr_limit != None else None
        if n == None and amp_limit == None and snr_limit == None:
            raise ValueError('No termination condition specified')

        i = 0
        while True:
            if n != None and i >= n:
                return i
            if self.count > 0:
                if amp_limit != None and self.amplitudes.min() < amp_limit:
                    return i
                if snr_limit != None and self.snr.min() < snr_limit:
                    return i
            self.step()
            i += 1
        return i

