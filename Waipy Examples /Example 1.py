import waipy
import numpy as np


z = np.linspace(0,2048,2048)
x = np.sin(50*np.pi*z)
y = np.cos(50*np.pi*z)

data_norm = waipy.normalize(x)
result = waipy.cwt(data_norm, 1, 1, 0.125, 2, 4/0.125, 0.72, 6, 
                   mother='Morlet',name='x')
waipy.wavelet_plot('Sine', z, data_norm, 0.03125, result)

data_norm1 = waipy.normalize(y)
result1 = waipy.cwt(data_norm1, 1, 1, 0.125, 2, 4/0.125, 0.72, 6, 
                    mother='Morlet',name='y')
waipy.wavelet_plot('Cosine', z, data_norm1, 0.03125, result1)

cross_power, coherence, phase_angle = waipy.cross_wavelet(result['wave'], 
                                                          result1['wave'])
waipy.plot_cross('Crosspower sine and cosine', cross_power, phase_angle, 
                 z, result, result1)
