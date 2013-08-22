from __future__ import absolute_import, division

import os, sys, pyfits
import numpy as np
from .utils import ExportDecorator

__author__ = 'Kolja Glogowski, Wiebke Herzberg'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'

__all__ = []
export = ExportDecorator(__all__)

def _auto_detect_format(fname, check_fits_type=True):
    """
    Auto-detects the file format of a given file.

    Auto-detection of the file format is done by first trying to open it as
    a FITS file. If this failes, an ascii file is assumed. If the file can
    be successfully opened as FITS file, the HDUs are (optionally) inspected
    to decide if the file contains a binary table or a plain image array.

    Parameters
    ----------
    fname : string
        File name.
    check_fits_type : bool (optional)
        Check for binary table FITS files.

    Returns
    -------
    fmt : string
        The file format: 'fits-img', 'fits-tbl' or 'ascii' if check_fits_type
        is True; otherwise 'fits' or 'ascii'.
    """
    try:
        hdus = pyfits.open(fname)
        if check_fits_type:
            fmt = 'fits-img'
            if len(hdus) > 1 and hdus[0].header.get('naxis') == 0 \
                             and hdus[1].header.get('xtension') == 'BINTABLE':
                fmt = 'fits-tbl'
        else:
            fmt = 'fits'
        hdus.close()
    except IOError:
        fmt = 'ascii'
    return fmt

@export
def read_timeseries(fname, fmt='auto'):
    """
    Reads time series from FITS or ASCII files.

    Parameters
    ----------
    fname : string
        The file name.
    fmt : string (optional)
        File format; possible values are 'auto', 'fits' or 'ascii'.

    Returns
    -------
    t : ndarray
        Time steps of the time series.
    a : ndarray
        Amplitude of the time series.
    """
    if fmt == 'auto':
        fmt = _auto_detect_format(fname)
        if fmt == 'fits-tbl':
            raise IOError('FITS tables are not supported for time series')
        elif fmt == 'fits-img':
            fmt = 'fits'
        assert fmt in [ 'fits', 'ascii' ]

    if fmt == 'fits':
        data = pyfits.getdata(fname).astype(np.float64)
    elif fmt == 'ascii':
        data = np.loadtxt(fname)
    else:
        raise ValueError('Invalid file format specified')

    if data.ndim != 2 or data.shape[1] != 2:
        raise IOError('Invalid number of columns')
    data = data[np.isfinite(data[:,1])]
    return data[:,0], data[:,1]

@export
def read_params(fname, fmt='auto', as_sarray=True, add_noise_cols=False):
    """
    Reads Pysca parameter files containing frequencies, amplitudes, phases,
    noise levels and signal-to-noise ratios.

    Parameters
    ----------
    fname : string
        The input file name.
    fmt : string (optional)
        File format; possible values are 'auto', 'fits' or 'ascii'.
    as_sarray : bool (optional)
        Set to False to return a plain ndarray instead of a structured array.
    add_noise_cols : bool (optional)
        Add noise/snr columns, even if the file does not contain noise/snr
        columns. Empty noise/snr columns are filled with nan values.

    Returns
    -------
    params : ndarray (plain or structured, depending on as_sarray)
        The parameters contained in the given file. If the result is a
        structured array, it contains the keys 'freq', 'amp', 'phase' and
        (if available) 'noise' and 'snr'. Missing columns are filled with
        nan values.
    """
    if fmt == 'auto':
        fmt = _auto_detect_format(fname, check_fits_type=False)
        assert fmt in [ 'fits', 'ascii' ]
    cnames = [ 'freq', 'amp', 'phase', 'noise', 'snr' ]
    if fmt == 'fits':
        data = pyfits.getdata(fname)
        if data.dtype.names == None:
            # Plain FITS file (ImageHDU).
            # Converting data to native float array with 3 or 5 columns.
            if data.ndim == 1:
                data = np.atleast_2d(data).T
            has_noise = (data.shape[1] >= 5)
            ncols = 5 if (has_noise or add_noise_cols) else 3
            a = np.empty((data.shape[0], ncols))
            a.fill(np.nan)
            n = min(ncols, data.shape[1])
            a[:,0:n] = data[:,0:n]
        else:
            # BinTable FITS file.
            # Creating plain array from bin table data with 3 or 5 columns.
            dnames = data.dtype.names
            if cnames[0] not in dnames:
                raise IOError("Missing 'freq' column")
            has_noise = (cnames[-1] in dnames) and (cnames[-2] in dnames)
            ncols = 5 if (has_noise or add_noise_cols) else 3
            a = np.empty((data.shape[0], ncols))
            a.fill(np.nan)
            for i in range(ncols):
                cni = cnames[i]
                if cni in dnames:
                    a[:,i] = data[cni]
    elif fmt == 'ascii':
        # ASCII file.
        # Adjusting array to 3 or 5 columns.
        data = np.atleast_2d(np.loadtxt(fname))
        has_noise = (data.shape[1] >= 5)
        ncols = 5 if (has_noise or add_noise_cols) else 3
        if data.shape[1] != 5:
            a = np.empty((data.shape[0], ncols))
            a.fill(np.nan)
            n = min(ncols, data.shape[1])
            a[:,0:n] = data[:,0:n]
        else:
            a = data
    else:
        raise ValueError('Invalid file format specified')

    if as_sarray:
        dnames = cnames[:ncols]
        sa = np.empty(a.shape[0], dtype=[(dni, 'f8') for dni in dnames])
        for i, dni in enumerate(dnames):
            sa[dni] = a[:,i]
        return sa
    else:
        return a

@export
def write_params(fname, params, fmt='fits-tbl', add_to_header=None,
                 clobber=False):
    """
    Writes Pysca parameter files containing frequencies, amplitudes, phases,
    noise levels and signal-to-noise ratios.

    Parameters
    ----------
    fname : string
        The output file name.
    params : ndarray (plain or structured)
        The params argument can be either a plain array with 3 or 5 columns,
        or a structured array with column names 'freq', 'amp', 'phase' and
        optionally 'noise' and 'snr'.
    fmt : string (optional)
        Supported formats are FITS binary tables ('fits-tbl'), (plain) FITS
        image files ('fits-img') and ASCII files ('ascii').
    add_to_header : None or list(pyfits.Card) (optional)
        Additional header entries. For ASCII files they are written as comment
        at the top of the file, right before the column header.
    clobber : bool (optional)
        Overwrite existing files.
    """
    params = np.asarray(params)
    if fmt == 'fits':
        fmt = 'fits-tbl'
    if fmt not in [ 'fits-tbl', 'fits-img', 'ascii' ]:
        raise ValueError('Invalid format string: %s' % fmt)

    # Check params array type.
    is_sarray = (params.dtype.names != None)
    if is_sarray:
        for cni in [ 'freq', 'amp', 'phase' ]:
            if cni not in params.dtype.names:
                raise ValueError('Missing column im params array: %s' % cni)
        has_noise = (('noise' in params.dtype.names) and
                     ('snr' in params.dtype.names))
    else:
        if params.ndim != 2:
            raise ValueError('Wrong params array dimension')
        if params.shape[1] != 3 and params.shape[1] != 5:
            raise ValueError('params array must have either 3 or 5 columns')
        has_noise = (params.shape[1] == 5)

    # Create output array.
    cnames = [ 'freq', 'amp', 'phase' ]
    if has_noise:
        cnames.extend([ 'noise', 'snr' ])
    if fmt == 'fits-tbl':
        a = np.zeros(params.shape[0], dtype=[(cni, 'f8') for cni in cnames])
        if is_sarray:
            for cni in cnames:
                a[cni] = params[cni]
        else:
            for i, cni in enumerate(cnames):
                a[cni] = params[:,i]
    else:  # fmt in [ 'fits-img', 'ascii' ]
        a = np.zeros((params.shape[0], len(cnames)), dtype='f8')
        if is_sarray:
            for i, cni in enumerate(cnames):
                a[:,i] = params[cni]
        else:
            for i, cni in enumerate(cnames):
                a[:,i] = params[:,i]

    # Check if the file exists and remove it if clobber is set.
    if os.path.exists(fname):
        if clobber:
            os.remove(fname)
        else:
            raise IOError("File '%s' already exists" % fname)

    # Write array to file.
    if fmt == 'ascii':
        # To be compatible with older versions of numpy.savetxt(), we write
        # the header manually and then pass the open file object to the
        # numpy.savetxt() routine.
        f = open(fname, 'wb' if sys.version_info[0] >= 3 else 'w')
        try:
            # Add additional header entries.
            if add_to_header:
                for ci in add_to_header:
                    f.write('# %s\n' % str(ci).strip())
            # Add column names to the header.
            f.write('# %s\n' % ' '.join(cnames).upper())

            # Use the numpy.savetxt() to write the array.
            np.savetxt(f, a)

        finally:
            f.close()

    else:  # fmt in [ 'fits-tbl', 'fits-img' ]:
        hdus = pyfits.HDUList()
        if fmt == 'fits-tbl':
            hdus.append(pyfits.BinTableHDU(a))
        else:  # fmt == 'fits-img'
            hdus.append(pyfits.PrimaryHDU(a))
            # Update NAXIS1 header entry; the way it is done should be
            # compatible with older pyfits versions.
            h = hdus[0].header
            h.update('naxis1', h['naxis1'], ', '.join(cnames).upper())

        # Add additional header entries. This could be done easier using the
        # header.append() routine, but this only works with PyFITS 3.x.
        if add_to_header != None:
            h = hdus[0].header
            for ci in add_to_header:
                if pyfits.__version__ < '3.1.2':
                    h.update(ci.key, ci.value, ci.comment)
                else:
                    h.set(ci.keyword, ci.value, ci.comment)

        # Finally write the FITS file.
        hdus.writeto(fname)

@export
def write_periodogram(fname, nu, p, add_to_header=None, clobber=False):
    """
    """
    nu, p = np.atleast_1d(nu, p)
    if nu.ndim > 1 or p.ndim > 1:
        raise ValueError('Input arrays must be 1d')
    if len(nu) != len(p):
        raise ValueError('Input arrays must have the same size')

    # Check if the file exists and remove it if clobber is set.
    if os.path.exists(fname):
        if clobber:
            os.remove(fname)
        else:
            raise IOError("File '%s' already exists" % fname)

    hdu = pyfits.PrimaryHDU(np.c_[nu, p])

    # Update NAXIS1 comment and add additional header entries. The way this
    # is done should be compatible with older PyFITS versions.
    h = hdu.header
    h.update('naxis1', h['naxis1'], 'NU, PERIODOGRAM')
    if add_to_header != None:
        for ci in add_to_header:
            if pyfits.__version__ < '3.1.2':
                h.update(ci.key, ci.value, ci.comment)
            else:
                h.set(ci.keyword, ci.value, ci.comment)

    # Write data to disk.
    hdu.writeto(fname)
