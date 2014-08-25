=====================
Wavelet.wavelet_plot
=====================
wavelet.wavelet_plot(time, data,dtmin, result)

        Open and read an archive .txt with data only (without date)

        :Parameters:    time:           array_like
                                        Raw of time sampling.

                        data:           number
                                        Time-sample of the vector. Example: Hourly, daily, monthly, etc...

			dtmin		number
					Power maximum resolution. Example: 0.03125

			result:         dict.

                                        All parameters for plot from wavelet.cwt.

			impath:		the path where you want to save the figures.		

        :Returns:       the plot

.. Seealso::
        wavelet.cwa_plot

Notes

This plot the wavelet transform!


Example::

	>> dt = 0.25

	>> date1 = 1871

	# Test data = sst_nino3.dat is already in the package!
	>> data,n,time = load_txt('sst_nino3.dat',dt,date1)

	# This normalize by variance
	>> data_norm, variance = normalize(data)

	# Continuous wavelet transform
	>> result = cwt(data_norm,0.25,variance,n,1,0.25,2*0.25,7/0.25,0.72,6,'Morlet')

	# Plot all results 
	>> wavelet_plot('SST_NINO3',time,data,0.03125,result,impath)

