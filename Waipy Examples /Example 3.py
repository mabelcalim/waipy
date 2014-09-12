import waipy
import numpy as np

# Create Signals
z = np.linspace(0,2048,2048)
x = np.sin(50 * np.pi * z) + 3.5 * np.random.randn(size(z))
y = np.cos(50 * np.pi * z + np.pi / 4) + 2.5 * np.random.randn(size(z))

#Norm and CWT for Sine with Noise:

data_norm = waipy.normalize(x)
result = waipy.cwt(data_norm, 1, 1, 0.25, 4, 4/0.25, 0.72, 6, 
                   mother='Morlet',name='x')
waipy.wavelet_plot('Sine with noise', z, data_norm, 0.03125, result)

# Norm and CWT for Cosine with Noise:
data_norm1 = waipy.normalize(y)
result1 = waipy.cwt(data_norm1, 1, 1, 0.25, 4, 4/0.25, 0.72, 6, 
                    mother='Morlet',name='y')
waipy.wavelet_plot('Cosine with noise', z, data_norm1, 0.03125, result1)

# Crosspower:
cross_power, coherence, phase_angle = waipy.cross_wavelet(result['wave'], 
                                                          result1['wave'])
waipy.plot_cross('Crosspower sine and cosine with Noise', cross_power, 
                 phase_angle, z, result, result1)
