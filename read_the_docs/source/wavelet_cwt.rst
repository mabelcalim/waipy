===================
Wavelet.cwt
===================
wavelet.cwt(data, dt, variance, n, pad, dj, s0, j1, lag1, param, mother)

        Continuous wavelet transform from data. Wavelet params can be modified as you wish.

        :Parameters:    data:           array_like.
                                        Raw of data or normalized data.
					
                        dt:             number.
                                        Time-sample of the vector. Example: Hourly, daily, monthly, etc...

                        variance:       number.
                                        Data variance.

			n:		number.
					Length of the data.

			pad:		number/flag.
					Pad the time series with zeroes to next pow of two length (recommended). 

					Default: pad = 1. 
			
			dj: 		number.
					Divide octave in sub-octaves. 
					If dj = 0.25 this will do 4 sub-octaves per octave.  

			s0:		number.
					The maximum frequency resolution. 
					If it is an annual data, s0 = 2*dt say start at a scale of 6 months.

					Default: s0 = 2*dt.

			j1:		number.
					Divide the power-of-teo with dj sub-octaves each.

				        Default: j1 =7/dj.	

			lag1:		number.
					Lag-1 autocoorelation for red noise background.

					Default: lag1 =0.72.	

			param: 		number/flag.
					The mother wavelet param.

					Default: param = 6 (Morlet function used as default).

			mother:		string.
					The mother wavelet funtion. 

					Default: moher = 'Morlet'.

        :Returns:       result:         dict.

                                        Return all parameters for plot.

.. Seealso::
        wavelet.cwa

Notes

The Morlet wavelet is used as default int this code. The wavelet.cwt function call all lib_wavelet.py functions:

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




Example ::

	>> dt = 0.25

	>> date1 = 1871

	# Test data = sst_nino3.dat is already in the package!

	>> data,n,time = load_txt('sst_nino3.dat',dt,date1)

	# This normalize by variance
	>> data_norm, variance = normalize(data)
	
	# Continuous wavelet transform
	>> result = cwt(data_norm,0.25,variance,n,1,0.25,2*0.25,7/0.25,0.72,6,'Morlet')

