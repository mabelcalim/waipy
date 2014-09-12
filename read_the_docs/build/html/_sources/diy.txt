How Waipy works
===============

This chapter shows on the Example: Niño3 SST how waipy calculates the CWT from a given time series.

CWT - Niño3 SST
---------------

0. Import python libraries
__________________________
::

    import numpy as np
    import pylab
    from pylab import detrend_mean
    import matplotlib.pyplot as plt
    import math


1. Load the data
________________

Before you can work with Waipy you have to load the data with eiter Waipy or another Python tool (eg. Pandas).

::

    # loading data with waipy.load_txt
    data,time = waipy.load_txt('sst_nino3.dat', 0.25, 1871)

.. automodule:: wavetest
    :members: load_txt


2. Normalize with standard score
________________________________


::

    # normalizing time series
    data_norm = waipy.normalize(data)

.. automodule:: wavetest
    :members: normalize

3. Choose Wavelet and parameters for the cwt
____________________________________________

::

    cwt(data, dt, pad, dj, s0, j1, lag1, param, mother, name):

.. automodule:: wavetest
    :members: cwt

The CWT calls now lib_wavelet.py where the calculations are made.

::

                                 +----------------+
                                 |    cwt.py 	  |
                                 +----------------+
                                        |
                                +----------------+
                                | lib_wavelet.py |
                                +----------------+
                                        |
                      +----------------+  +----------------+
                      |  def wavelet   |--| def wave_signif|
                      +----------------+  +----------------+
                              |
            +----------------+  +----------------+
            | def nextpow2   |--| def wave_bases |
            +----------------+  +----------------+

4. Calculate the daughter wavelets:
___________________________________
::

    # lib_wavelet.wavelet(data, dt, param, dj, s0, j1, mother)

.. automodule:: lib_wavelet
    :members: wavelet

::

        n1 = len(Y)  # time series length 
        s0 = 2*dt  # smallest scale of the wavelet
        dj = 0.25  # spacing between discrete scales


- J1 :doc:`checkitout`

::

    J1= int(np.floor((np.log10(n1*dt/s0))/np.log10(2)/dj))  # J1+1 total os scales


- Call nextpow2 :doc:`checkitout`

::

    if (pad ==1):
        base2 = nextpow2(n1) #call det nextpow2
        n = base2

- k :doc:`checkitout`::

        # construct wavenumber array used in transform
        # simetric eqn 5  
        k = np.arange(n/2)
        import math
        k_pos,k_neg=[],[]
        for i in range(0,n/2+1):
                k_pos.append(i*((2*math.pi)/(n*dt)))  # frequencies as in eqn5
                k_neg = k_pos[::-1]  # inversion vector
                k_neg = [e * (-1) for e in k_neg]                       # negative part 
                k_neg = k_neg[1:-1]  # delete the first value of k_neg = last value of k_pos
        k = np.concatenate((k_pos,k_neg), axis =1)  # vector of symmetric


::

        +---------------------------------------+
        |                 1                     |
        |             +------+                  | 
        |             | 0.0  |                  |
        |             |      |                  |
        |len(k) = 512 |      |                  |
        |             |      |                  |
        |             |-0.049|                  |
        |             +------+                  |
        +---------------------------------------+

- f :doc:`checkitout`::
                         
        # compute fft of the padded time series
        f = np.fft.fft(x,n)
       
::

        +---------------------------------------+
        |                 1                     |
        |             +------------+            | 
        |             |-9.33e-15+0j|            |
        |             |            |            |
        |len(k) = 512 |            |            |
        |             |            |            |
        |             |31.9 -5.55j |            |
        |             +------------+            |
        +---------------------------------------+
      

5. Wavebases:
_______________________________
::

    # daughter, fourier_factor, coi, dofmin = wave_bases(mother, k, scale[a1], param)

.. automodule:: lib_wavelet
    :members: wave_bases

- expnt :doc:`checkitout`::

    expnt = -pow(scale*k-k0,2)/2*(k>0) 

.. note::
    Only the values of the last scale are avaiable on checkitout.



::

    +---------------------------------------+
    |               scale[a1] = 32          |
    |             +-----------------+       |
    |             |-0               |       |
    |             |                 |       |
    |len(k) = 512 |                 |       |
    |             |                 |       |
    |             |               -0|       |
    |             +-----------------+       |
    +---------------------------------------+

- norm :doc:`checkitout`::

    norm = math.sqrt(scale*k[1])*(pow(math.pi,-0.25))*math.sqrt(len(k))

::

    +---------------------------------------+
    |          scale[a1] = 32               |
    |       +-----------------+             |
    |       |2.66        39.07|             |
    |       +-----------------+             |
    |                                       |
    +---------------------------------------+

- daughter :doc:`checkitout`::

                daughter = []  # define daughter as a list
                for ex in expnt:  # for each value scale (equal to next pow of 2)
                        daughter.append(norm*math.exp(ex))
                k = np.array(k)  # turn k to array
                daughter = np.array(daughter)  # transform in array
                daughter = daughter*(k>0)  # Heaviside step function

.. note:: Only the values of the last scale are avaiable on checkitout.

::

    +---------------------------------------+
    |               scale[a1] = 32          |
    |             +-----------------+       | 
    |             |0.0              |       |
    |             |                 |       |
    |len(k) = 512 |                 |       |
    |             |                 |       |
    |             |              0.0|       |
    |             +-----------------+       |
    +---------------------------------------+

6. Significance of Wavelet:
_____________________________
::

    # signif, fft_theor = lib_wavelet.wave_signif(1.0, dt, scale, 0, lag1, 0.95, -1, mother, param)

.. automodule:: lib_wavelet
    :members: wave_signif

7. FFT of data and Levels:
______________________________________

.. automodule:: wavetest
    :members: fft

.. automodule:: wavetest
    :members: levels


8. CWT Plot:
__________________

.. automodule:: wavetest
    :members: wavelet_plot

