Installation
============



**1) Download the package --> waipy-version.tar.gz**

https://pypi.python.org/pypi/waipy/

**2) Unpack it in you download directory**

::

    tar -vzxf waipy-version.tar.gz

**3) Enter the new directory**

:: 

    cd waipy-version

**4) Install package**

::

    python setup.py install

or::

    pip install -e .

**5) Open ipython or IpythonNotebook and copy the following for a test:**

::

    import waipy
    # loading data for test
    data,time = waipy.load_txt('sst_nino3.dat', 0.25, 1871)
    # normalizing time series
    data_norm = waipy.normalize(data)
    # calculating continuos wavelet transform using Morlet wavelet
    result = waipy.cwt(data_norm, 0.25, 1, 0.25, 2*0.25, 7/0.25, 0.72, 6, 'Morlet')
    # ploting the result - impath: 'Users/you/figures/'
    waipy.wavelet_plot('SST_NINO3', time, data, 0.03125, result, impath)


    # CROSS ANALYSIS
    cross_power, coherence = waipy.cross_wavelet(result['wave'], result['wave'])
    waipy.plot_cross(cross_power, time, result)
    waipy.plot_cohere(coherence, time, result)


