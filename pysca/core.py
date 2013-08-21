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
        self._snr_width = float(snr_width) if snr_width else None
        self._ofac = float(ofac)
        self._hifreq = float(hifreq) if hifreq != None else self._numax
        self._verbose = int(verbose)
        self._vprint = VerbosePrinter(self._verbose)
        self._prev_ts = self._next_ts = None
        self._nu = self._orig_per = self._prev_per = self._next_per = None

        dnames = [ 'freq', 'amp', 'phase', 'noise', 'snr' ]
        if not self._snr_width:
            dnames = dnames[:3]
        self._params = np.empty(0, dtype=[(dni, 'f8') for dni in dnames])
        self._params = self._params.view(np.recarray)

    def _append_to_params(self, freq, amp, phase, noise=np.nan, snr=np.nan):
        a = np.empty(self._params.shape[0]+1, dtype=self._params.dtype)
        a[:-1] = self._params
        if len(a.dtype) == 3:
            a[-1] = (freq, amp, phase)
        else:
            a[-1] = (freq, amp, phase, noise, snr)
        self._params = a.view(np.recarray)

    def _update_params(self, freq, amp, phase, noise=None, snr=None):
        a = np.empty(self._params.shape[0]+1, dtype=self._params.dtype)
        a = a.view(np.recarray)
        a.freq = freq
        a.amp = amp
        a.phase = phase
        if len(a.dtype) > 3:
            if noise != None:
                a.noise = noise
            else:
                a.noise.fill(np.nan)
            if snr != None:
                a.snr = snr
            else:
                a.snr.fill(np.nan)
        self._params = a

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
        return len(self._params)

    @property
    def result(self):
        return self._params

    def _calc_periodogram(self, t, a):
        return utils.compute_periodogram(t, a, self.ofac, self.hifreq)

    def _find_highest_peak(self, nu, per, use_nuidx=True):
        if use_nuidx:
            # Limit the periodogram to selected frequency range.
            nu = nu[self._nuidx]
            per = per[self._nuidx]
        return utils.find_highest_peak(nu, per)

    def step(self):
        vprint = self._vprint
        if self._next_per == None:
            self._next_per = self.orig_periodogram
            self._next_ts = self.orig_ts

        # Select the next periodogram
        nu, per = self.nu, self.next_periodogram

        # Find frequency of the highest peak in the current periodogram
        new_freq = np.r_[self._params.freq, self._find_highest_peak(nu, per)]

        # Fit original time series using already extracted mode parameters
        new_amp, new_phase, ok, misc = fit.fit_timeseries(self.t,
            self.orig_ts, new_freq, self._params.amp, self._params.phase)
        if not ok:
            raise PyscaError('Harmonic fit failed for peak frequency %f' % (
                new_freq[-1]))

        # Prewhiten the original time series using the new mode parameters
        new_ts = fit.prewhiten(self.t, self.orig_ts, new_freq, new_amp,
                               new_phase)

        # Compute new periodogram from the prewhitened time series
        new_nu, new_per = self._calc_periodogram(self.t, new_ts)

        # Compute noise for the last extracted frequency using the median
        if self.snr_width:
            noise = utils.median_noise_level(
                nu, new_per, new_freq[-1], self.snr_width)
            snr = new_amp[-1] / noise
            new_noise = np.r_[self._params.noise, noise]
            new_snr = np.r_[self._params.snr, snr]
        else:
            new_noise = new_snr = None

        if self.snr_width:
            vprint('freq = %f, amp = %f, phase = %f, noise = %f, snr = %f' % (
                    new_freq[-1], new_amp[-1], new_phase[-1],
                    new_noise[-1], new_snr[-1]), v=1)
        else:
            vprint('freq = %f, amp = %f, phase = %f' % (
                    new_freq[-1], new_amp[-1], new_phase[-1]), v=1)

        # Update object data
        self._update_params(new_freq, new_amp, new_phase, new_noise, new_snr)
        self._prev_per = self._next_per
        self._prev_ts = self._next_ts
        self._next_per = new_per
        self._next_ts = new_ts

    def check_term_conditions(self, amp, snr):
        """
        Returns
        -------
        True, if the limits are met; otherwise False.
        """
        if len(self._params) == 0:
            return True
        amp = float(amp) if amp != None else None
        snr = float(snr) if snr != None else None
        if amp != None and self._params.amp[-1] < amp:
            return False
        if snr != None and self._params.snr[-1] < snr:
            return False
        return True

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

        Returns
        -------
            number of iterations
        """
        n = int(n) if n != None else None
        if n == None and amp_limit == None and snr_limit == None:
            raise ValueError('No termination condition specified')

        i = 0
        while True:
            if n != None and i >= n:
                return i
            if not self.check_term_conditions(amp_limit, snr_limit):
                return i
            self.step()
            i += 1
        return i
