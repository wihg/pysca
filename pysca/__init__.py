from __future__ import absolute_import

from .core import *
from .fit import *
from .utils import *
from . import core, fit, utils

__all__ = []
__all__.extend(core.__all__)
__all__.extend(fit.__all__)
__all__.extend(utils.__all__)

__author__ = 'Wiebke Herzberg, Kolja Glogowski'
__maintainer__ = 'Kolja Glogowski'
__email__ = 'kolja@kis.uni-freiburg.de'
__license__ = 'MIT'
__version__ = '0.0.2'
