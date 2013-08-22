from pysca import Pysca
from pysca.io import read_timeseries, write_params

t, a = read_timeseries('timeseries.fits')

p = Pysca(t, a, 0.5, 20, 3, verbose=1)
p.run(5)

write_params('out.dat', p.result, fmt='ascii', clobber=True)

