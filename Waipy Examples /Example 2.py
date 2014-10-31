import waipy
import numpy as np

x = np.linspace(0, 100, 100)
y1 = np.random.rand(100)  # Generation of the Random Signal 1
y2 = np.random.rand(100)  # Generation of the Random Signal 2

data_norm = waipy.normalize(y1)
data_norm1 = waipy.normalize(y2)

result = waipy.cwt(data_norm, 1, 1, 0.25, 2, 4/0.25, 0.72, 6, 
                   mother='Morlet', name='x')
result1 = waipy.cwt(data_norm1, 1, 1, 0.25, 2, 4/0.25, 0.72, 6, 
                    mother='Morlet', name='y')

waipy.wavelet_plot('y1-random-signal', x, data_norm, 0.03125, result)
waipy.wavelet_plot('y2-random-signal', x, data_norm1, 0.03125, result1)

cross_power, coherence, phase_angle = waipy.cross_wavelet(result['wave'], 
                                                          result1['wave'])
waipy.plot_cross('Crosspower of y1 and y2', cross_power, phase_angle, x, 
                 result, result1)

