DIY 
===

Do it Yourself (DIY) your own wavelet toolkit!

CWT - Niño3 SST
---------------

0. Import python libraries
__________________________
::

	import numpy as np
	import pylab
	from pylab import detrend_mean
	import math


1.Choose and implement the wavelet function
____________________________________________

    
 ::     

        
	def wave_bases(mother,k,scale,param):
          """Computes the wavelet function as a function of Fourier frequency
           used for the CWT in Fourier space (Torrence and Compo, 1998) 
           -- This def is called automatically by def wavelet --

           _____________________________________________________________________
           Inputs:
                  mother - a string equal to 'Morlet'
                  k      - a vectorm the Fourier frequecies 
                  scale  - a number, the wavelet scale
                  param  - the nondimensional parameter for the wavelet function

           Outputs:
                  daughter       - a vector, the wavelet function 
                  fourier_factor - the ratio os Fourier period to scale
                  coi            - a number, the cone-of-influence size at the scale
                  dofmin         - a number, degrees of freedom for each point in the wavelet power (Morlet = 2)
           
           Call function:
                  daughter,fourier_factor,coi,dofmin = wave_bases(mother,k,scale,param) 
           _____________________________________________________________________

.. note::
    The Morlet wavelet is used as default in this code.


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

                daughter = []                                           # define daughter as a list
                for ex in expnt:                                        # for each value scale (equal to next pow of 2)
                        daughter.append(norm*math.exp(ex))
                k = np.array(k)                                         # turn k to array
                daughter = np.array(daughter)                           # transform in array
                daughter = daughter*(k>0)                               # Heaviside step function
.. note::
    Only the values of the last scale are avaiable on checkitout.    




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


2. Find the next pow of 2 
___________________________
 ::

  def nextpow2(i):
      n = 2
      while n < i: n = n * 2
      return n



3. Compute wavelet power spectrum
_________________________________
 ::

  def wavelet(Y,dt,mother,param):#,pad,dj,s0,J1,mother,param):
        """Computes the wavelet continuous transform of the vector Y, by definition:
        
           W(a,b) = sum(f(t)*psi[a,b](t) dt)               a dilate/contract 
           psi[a,b](t) = 1/sqrt(a) psi(t-b/a)              b displace 

           Only Morlet wavelet (k0=6) is used
           The wavelet basis is normalized to have total energy = 1 at all scales 
           _____________________________________________________________________
           Input:
                Y - time series 
                dt - sampling rate
                mother - the mother wavelet function
                param - the mother wavelet parameter

           Output:
                ondaleta - wavelet bases at scale 10 dt
                wave - wavelet transform of Y
                period - the vector of "Fourier"periods ( in time units) that correspond to the scales
                scale - the vector of scale indices, given by S0*2(j*DJ), j =0 ...J1
                coi - cone of influence                   
      
           Call function:
           ondaleta, wave, period, scale, coi = wavelet(Y,dt,mother,param)  
           _____________________________________________________________________

::

        n1 = len(Y)                                                     # time series length 
        s0 = 2*dt                                                       # smallest scale of the wavelet
        dj = 0.25                                                       # spacing between discrete scales


- J1 :doc:`checkitout`
::

	J1= int(np.floor((np.log10(n1*dt/s0))/np.log10(2)/dj))          # J1+1 total os scales


- Call nextpow2 :doc:`checkitout`
::

        if (pad ==1) :
                base2 = nextpow2(n1)                                    #call det nextpow2
        n = base2

- k :doc:`checkitout`::

        # construct wavenumber array used in transform
        # simetric eqn 5  
        k = np.arange(n/2)
        import math
        k_pos,k_neg=[],[]
        for i in range(0,n/2+1):
                k_pos.append(i*((2*math.pi)/(n*dt)))                    # frequencies as in eqn5
                k_neg = k_pos[::-1]                                     # inversion vector
                k_neg = [e * (-1) for e in k_neg]                       # negative part 
                k_neg = k_neg[1:-1]                                     # delete the first value of k_neg = last value of k_pos
        k = np.concatenate((k_pos,k_neg), axis =1)                      # vector of symmetric


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
      
