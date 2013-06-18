from __future__ import absolute_import

import collections
import numpy as np
from ._fasper import fasper
from ._prewhiten import prewhiten

__author__ = 'Wiebke Herzberg, Kolja Glogowski'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

class Pysca(object):
    def __init__(self, t, a, numin, numax, snr_width, ofac=8.0, hifreq=None):
        """
        Pysca(t, a, numin, numax, snr_width, ofac=8.0, hifreq=None)

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

        """
        self._t = np.asarray(t, dtype=np.float64)
        self._a = np.asarray(a, dtype=np.float64)
        if self._t.ndim != 1 or self._a.ndim != 1:
            raise ValueError('Input arrays must be 1d.')
        elif self._t.size != self._a.size:
            raise ValueError('Input arrays must have the same sizes.')
        self._numin, self._numax = float(numin), float(numax)
        self._snr_width = float(snr_width) if snr_width != None else None
        self._ofac = float(ofac)
        self._hifreq = float(hifreq) if hifreq != None else self._numax
        self._last_ts = self._next_ts = None
        self._nu = self._orig_per = self._last_per = self._next_per = None
        self._freqs = []     # List of extracted frequencies
        self._noise = []     # List of corresponding noise levels
        self._optpar = None  # Start values (am, ph) for the fit (prewhiten)

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
    def last_ts(self):
        return self._last_ts

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
    def last_periodogram(self):
        return self._last_per

    @property
    def next_periodogram(self):
        return self._next_per

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
        nu, p, nout, jmax, var = fasper(t, a, self.ofac, hifreq=self.hifreq)
        # Note: amplitudes from fasper are not scaled correctly. To get
        # amplitudes corresponding to the amplitude values from the fit, do
        # the normalization 2.0*sqrt(var*wk2(i)/n) (from kperiod).
        psq = 2.0 * np.sqrt(var * p / len(a))
        return nu, psq

    def _find_highest_peak(self, nu, per, use_nuidx=True):
        if use_nuidx:
            # Limit the periodogram to selected frequency range
            nu = nu[self._nuidx]
            per = per[self._nuidx]

        # Get index of highest peak in selected range
        imax = np.argmax(per)

        # Determine the frequency value by parabolic interpolation, which is
        # only possible if the maximum is not at the beginning of or end
        # of the periodogram.
        if imax == 0 or imax == per.size - 1:
            peak = per[imax]
        else:
            # Get values around the maximum
            frq1 = nu[imax-1]
            frq2 = nu[imax]
            frq3 = nu[imax+1]
            y1 = per[imax-1]
            y2 = per[imax]
            y3 = per[imax+1]

            # Parabolic interpolation formula
            t1 = (y2-y3) * (frq2-frq1)**2 - (y2-y1) * (frq2-frq3)**2
            t2 = (y2-y3) * (frq2-frq1) - (y2-y1) * (frq2-frq3)
            peak = frq2 - 0.5 * t1/t2
        return peak

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
        new_ts, new_optpar = prewhiten(self.t, self.orig_ts, self._freqs,
                                       self._optpar)

        # amplitude of the last extracted frequency
        amp = new_optpar[-1][0]

        # Compute new periodogram from the prewhitened time series
        new_nu, new_per = self._calc_periodogram(self.t, new_ts)

        # Calculate noise for the last extracted frequency using the median
        idx = (nu >= freq - 0.5 * self.snr_width) & (
               nu <= freq + 0.5 * self.snr_width)
        noise = np.median(new_per[idx])
        self._noise.append(noise)

        # Update object data
        self._last_per = self._next_per
        self._last_ts = self._next_ts
        self._next_per = new_per
        self._next_ts = new_ts
        self._optpar = new_optpar


    def run(self, n=None, amp_limit=None, snr_limit=None):
        """
        Parameters
        ----------
        n : int or None
            number of steps; default: no limit
        amp_limit : float or None
            amplitude limit (termination condition); default: no limit
        snr_limit : float or None
            signal-to-noise limit (termination condition); default: no limit

        """
        n = int(n) if n != None else None
        amp_limit = float(amp_limit) if amp_limit != None else None
        snr_limit = float(snr_limit) if snr_limit != None else None

