A tool for automated frequency extraction from photometric time series.

Pysca is a software tool that allows automated extraction of frequencies, amplitudes and phases from non-equally sampled time series. The extraction is done by identifying the highest peaks in the Lomb-Scargle periodogram and fitting the time series with a sum of harmonic functions of the corresponding frequencies. This is implemented using an iterative algorithm where the time series is progressively prewhitened up to an user defined termination condition. The signal-to-noise ratio is calculated for every frequency as a statistical measure of significance.