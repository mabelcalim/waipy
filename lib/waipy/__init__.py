# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 2013

@author: Mabel Calim Costa
"""

from __future__ import absolute_import

from .cwt.lib_wavelet import nextpow2
from .cwt.lib_wavelet import wave_bases
from .cwt.lib_wavelet import wavelet
from .cwt.lib_wavelet import wave_signif
from .cwt.wavetest    import load_txt
from .cwt.wavetest    import load_nc
from .cwt.wavetest    import normalize
from .cwt.wavetest    import cwt
from .cwt.wavetest    import wavelet_plot
from .cwt.wavetest    import fft

#from .cwt import wavetest


from .cwa.cross_wavelet import cross_wavelet
from .cwa.cross_wavelet import plot_cross
from .cwa.cross_wavelet import plot_cohere



# Define the objects imported by imports of the form: from pyclimatetools import *
__all__ = ['nextpow2', 'wave_bases', 'wavelet', 'wave_signif','load_txt','load_nc','normalize','cwt','wavelet_plot','fft','cross_wavelet','plot_cross', 'plot_cohere']

# Package version number.
__version__ = '0.0.9.0'
