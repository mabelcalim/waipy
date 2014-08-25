===================
Wavelet.load_nc
===================
wavelet.load_txt(file,var,dt,date1)

        Open and read an archive .txt with data only (without date)

        :Parameters:    file:	        vector
                                        Vector archive to be open.
			var:		string
					variable from archive.nc
                        dt:             number
                                        Time-sample of the vector. Example: Hourly, daily, monthly, etc...
                        date1:          number
                                        The initial time of the data. Example: 1985.

        :Returns:       data:           array_like
                                        Raw of data
                        time:           array_like
                                        Raw of time sampling.

.. Seealso::
        wavelet.load_txt




Example::
	
	# Creating a netcdf file 

	# download the pupynere package

	# more info : https://bitbucket.org/robertodealmeida/pupynere/

	>> import pupynere as pp

 	>> f = pp.netcdf_file('simple.nc', 'w')
	
	>> f.history = 'Created for a test'

	>> f.createDimension('time', 10)

	>> time = f.createVariable('time', 'i', ('time',))

	>> time[:] = range(10)

	>> time.units = 'days since 2008-01-01'

	>> f.close()

        # load netcdf

	>> dt = 1

	>> data, time = waipy.load_nc('simple.nc','time',dt,2000)


