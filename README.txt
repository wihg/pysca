=====
Pysca
=====

Pysca is a software package that allows automated extraction of frequencies, amplitudes and phases from non-equally sampled time series of heat-driven pulsators. The extraction is carried out by identifying the highest peaks in the Lomb-Scargle periodogram and fitting the time series with a superposition of harmonic functions of the corresponding frequencies. This is implemented using an iterative algorithm where the time series is progressively prewhitened up to an user defined termination condition. The signal-to-noise ratio is calculated for every frequency as a statistical measure of significance.

Requirements
============

Pysca runs with Python version 2.6 and 2.7; there is currently no support for Python 3. It depends on the following Python packages:

- NumPy, version 1.4.1 or newer
- SciPy, version 0.7.2 or newer
- PyFITS, version 2.3.1 or newer

If you use the ``pip`` package manager (see below), these packages will be (if neccessary) automatically downloaded and installed. Note however that ``pip`` will install these packages from source and that the build process might take quite some time (especially for the rather large SciPy package). So it might be a good idea to install a pre-compiled version of these packages (e.g. by using the package manager of your Linux distribution) before you install Pysca.

Installation
============

The prefered way to install Pysca is to use the `pip <https://pypi.python.org/pypi/pip>`_::

    pip install --user pysca

This way Pysca will be installed completely into the user's home directory. The ``pysca`` executable will be installed to::

      $HOME/.local/bin/pysca                       (Linux, MacOS w/ Python 2.6)
  or  $HOME/Library/Python/X.Y/bin/pysca                  (MacOS w/ Python 2.7)

where ``X.Y`` is the Python version (i.e. ``2.7``) and the ``pysca`` library will be
installed under::

      $HOME/.local/lib/pythonX.Y/site-packages     (Linux, MacOS w/ Python 2.6)
  or  $HOME/Library/Python/X.Y/lib/python/site-packages   (MacOS w/ Python 2.7)

Another advantage of ``pip`` is, that you can easily upgrade Pysca to the latest version, by entering::

    pip install --user -U --no-deps pysca

or completely uninstall Pysca using the following command::

    pip uninstall pysca

To run the ``pysca`` executable, you have to make sure, that it is in your ``PATH``, by either adding the directory, where the executable was installed, to your ``PATH`` environment variable, e.g. using the *Bash*::

    export PATH=$HOME/.local/bin/pysca:$PATH

or by creating a symlink of the binary to a directory which is already included in your ``PATH`` environment variable, e.g.::

    ln -s $HOME/.local/bin/pysca $HOME/bin/

If you, for any reason, don't want to use ``pip`` or any other Python installer (like `easy_install`), you can also download the ``.tar.gz`` file, extract it to some directory and then source the provided ``env_pysca.sh`` script::

    tar -xzf pysca-x.y.z.tar.gz
    source pysca-x.y.z/env_pysca.sh

Note however, that this script only works with the *Bash*; if you are using the *C shell* or anyting else, you can still use Pysca without using ``pip``, but you have to set the ``PYTHONPATH`` and ``PATH`` environment variables by your self.

How to use
==========

To use Pysca, you need an input file containing the (detrended) time series, either in FITS or in ASCII format. The file must have two columns, with the time steps in the first one and the amplitudes in the second one. Note that the physical unit of the time column determines the unit of the resulting frequencies. If the time column contains fractional days (e.g. JD or MJD), then the frequency unit is c/d (cycles per day). If, on the other hand, the time column contains values in seconds, then all frequencies are in Hz.

The ``pysca`` command
---------------------
If Pysca was installed correctly, you can start it with the ``pysca`` command. An overview of all command-line options is available using the ``pysca --help`` command. In order to extract mode parameters from the time series, it is required that you specify a frequency range, which restricts the part of the spectrum where mode parameters should be extracted from (option ``-f``). In addition you need to use the ``-b`` option to specify the width around each peak that is used to estimate the SNR; you can disable the SNR computation by explicitely setting it to 0. Finally you have to specify the output file name using the ``-o`` option.

The following gives a minimal example on how to use the ``pysca`` command::

    pysca timeseries.fits -o out.dat -f 0.5 20 -b 3 -n 5

This reads the time series from a file called ``timeseries.fits``, extracts the mode parameters of the 5 highest peaks from the frequency range between 0.5 and 20 c/d in the periodogram (assuming the time column of the time series contains MJD values), uses a width of 3 c/d around each peak to estimated the noise level and finally writes the result to an ASCII file called ``out.dat``.


The python interface
--------------------
As an alternative to the ``pysca`` command, you can also use the Python interface of Pysca. The following Python example does essentially the same as the example above::

    from pysca import Pysca
    from pysca.io import read_timeseries, write_params

    t, a = read_timeseries('timeseries.fits')

    p = Pysca(t, a, 0.5, 20, 3, verbose=1)
    p.run(5)

    write_params('out.dat', p.result, fmt='ascii', clobber=True)

**Important note:** Pysca is still in an early development state, so it is most likely that some of the command-line options, aswell as parts of the Python API will change in the future.
