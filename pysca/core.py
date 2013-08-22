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
                 params=None, verbose=0):
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
        hifreq : None or float
            maximum frequency of the Lomb-Scargle periodogram; default: numax
        params : None or structured array
            custom parameters, e.g. from a previous run, default: None
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
        if params != None and (dnames[0] in params.dtype.names):
            self._params = np.empty(
                len(params), dtype=[(dni, 'f8') for dni in dnames])
            for dni in dnames:
                if dni in params.dtype.names:
                    self._params[dni] = params[dni]
                else:
                    self._params[dni] = np.nan
        else:
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
        self._vprint('Computing periodogram...', v=2)
        return utils.compute_periodogram(t, a, self.ofac, self.hifreq)

    def _find_highest_peak(self, nu, per, use_nuidx=True):
        if use_nuidx:
            # Limit the periodogram to selected frequency range.
            nu = nu[self._nuidx]
            per = per[self._nuidx]
        return utils.find_highest_peak(nu, per)

    def step(self):
        vprint = self._vprint

        # Check if the needed periodogram is in place
        if self._next_per == None:
            if len(self._params) == 0:
                # No parameters yet, start with the original periodogram
                vprint('Initializing...', v=1)
                self._next_per = self.orig_periodogram
                self._next_ts = self.orig_ts
            else:
                # There are already parameters, so we need to prewhiten the
                # original time series, in order to get the next periodogram;
                # in case any amplitudes or phases are missing or wrong, we
                # first perform a fit to the time series.
                vprint('Initializing [nfreqs: %d]...' % len(self._params), v=1)
                vprint('Fitting time series...', v=2)
                self._params.amp[np.isnan(self._params.amp)] = 1.0
                self._params.phase[np.isnan(self._params.phase)] = 0.5
                new_amp, new_phase, ok, misc = fit.fit_timeseries(self.t,
                    self.orig_ts, self._params.freq, self._params.amp,
                    self._params.phase)
                if not ok:
                    raise PyscaError('Harmonic fit failed [ier=%d]' % misc[3])
                vprint('Prewhitening...', v=2)
                new_ts = fit.prewhiten(self.t, self.orig_ts, self._params.freq,
                                       new_amp, new_phase)
                new_nu, new_per = self._calc_periodogram(self.t, new_ts)
                self._next_per = new_per
                self._next_ts = new_ts
                self.orig_periodogram
            vprint(v=1)

        # Select the next periodogram
        nu, per = self.nu, self.next_periodogram

        # Find frequency of the highest peak in the current periodogram
        new_freq = np.r_[self._params.freq, self._find_highest_peak(nu, per)]
        vprint('Frequency #%d:' % len(new_freq), new_freq[-1], v=1)

        # Fit original time series using already extracted mode parameters
        vprint('Fitting time series...', v=2)
        new_amp, new_phase, ok, misc = fit.fit_timeseries(self.t,
            self.orig_ts, new_freq, self._params.amp, self._params.phase)
        if not ok:
            raise PyscaError('Harmonic fit failed for peak frequency %f' % (
                new_freq[-1]) + ' [ier=%d]' % misc[3])
        vprint('Freq: ', new_freq[-1], ', Amp: ', new_amp[-1], ', Phase: ',
               new_phase[-1], sep='', v=1)

        # Prewhiten the original time series using the new mode parameters
        vprint('Prewhitening...', v=2)
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
            vprint('Noise: ', new_noise[-1], ', SNR: ',new_snr[-1],
                   sep='', v=1)
        else:
            new_noise = new_snr = None

        # Update object data
        self._update_params(new_freq, new_amp, new_phase, new_noise, new_snr)
        self._prev_per = self._next_per
        self._prev_ts = self._next_ts
        self._next_per = new_per
        self._next_ts = new_ts

    def check_term_conditions(self, amp_limit=None, snr_limit=None):
        """
        Check if the last extracted mode parameters are above the provided
        limits.

        Parameters
        ----------
        amp_limit : scalar
            amplitude limit; default: no limit
        snr_limit : scalar
            signal-to-noise limit; default: no limit

        Returns
        -------
        True, if the limits are met; otherwise False.
        """
        if len(self._params) == 0:
            return True
        amp_limit = float(amp_limit) if amp_limit != None else None
        snr_limit = float(snr_limit) if snr_limit != None else None
        if amp_limit != None and self._params.amp[-1] < amp_limit:
            return False
        if snr_limit != None and self._params.snr[-1] < snr_limit:
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
        vprint = self._vprint
        n = int(n) if n != None else None
        if n == None and amp_limit == None and snr_limit == None:
            raise ValueError('No termination condition specified')

        i = 0
        while True:
            if n != None and i >= n:
                break
            if not self.check_term_conditions(amp_limit, snr_limit):
                break
            self.step()
            vprint(v=1)
            i += 1
        return i
