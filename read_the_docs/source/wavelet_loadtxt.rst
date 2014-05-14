===================
Wavelet.load_txt
===================
wavelet.load_txt(archive,dt,date1)

	Open and read an archive .txt with data only (without date)

	:Parameters: 	archive : 	vector
				 	Vector archive to be open.
			dt:		number
				 	Time-sample of the vector. Example: Hourly, daily, monthly, etc...
			date1:		number
					The initial time of the data. Example: 1985.

	:Returns:	data:		array_like
				 	Raw of data	
			n:		number
				        The length of the data
			time:		array_like
					Raw of time sampling.

.. Seealso::
	wavelet.load_netcdf

Notes


This function is linked to data/txt directory, so, if you have any file extension .txt put it in the following folder:
::
/lib/wavelet/data/txt


Example::

	>> dt = 0.25

	>> date1 = 1871

	# Test data = sst_nino3.dat is already in the package!

	>> data,n,time = load_txt('sst_nino3.dat',dt,date1)

		
					
