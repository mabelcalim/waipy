===================
Wavelet.normalize
===================
wavelet.normalize(data)

        Data normalize by variance. The mean value is removed. 

        :Parameters:    data:           array_like
                                        Raw of data

        :Returns:       data:           array_like
                                        Raw of data normalized

                        variance:       Data variance
 

Notes

You can skip this function if it the normalization is not necessary (e.g. EOF data).


Example::

	>> dt = 0.25

	>> date1 = 1871

	# Test data = sst_nino3.dat is already in the package!

	>> data,n,time = load_txt('sst_nino3.dat',dt,date1) 

	# This normalize by variance
	>> data_norm, variance = normalize(data)
