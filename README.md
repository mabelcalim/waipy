waipy
=====
This guide includes a Continuous Wavelet Transform (CWT), significance  tests
from based on Torrence and Compo (1998)
![ScreenShot](https://wavelet-analysis.readthedocs.org/en/latest/_images/nino_wavelet.png)


and Cross Wavelet Analysis  (CWA) based on Maraun and Kurths(2004).
![ScreenShot](https://wavelet-analysis.readthedocs.org/en/latest/_images/salt_OGCM_cross.png)

Installation
============

A local install can be done using the provided setup.py file:

    python setup.py install

The installation path can be changed using the **--prefix** switch, e.g.:

    python setup.py install --prefix $HOME/inst

If you plan on modifying the code, use the **develop** target in combination
with the **--user** swich:

    export PYTHONUSERBASE=$HOME/inst/pip_installs
    export PYTHONPATH=/home/mweigand/inst/pip_installs/lib/python2.7/\
    site-packages/:$PYTHONPATH
    export PATH=$HOME/inst/pip_installs/bin:$PATH
    python setup.py develop --user

The first three lines should also be included in the **$HOME/.bashrc** file.

Requirements
============

* netcdf4 (https://github.com/Unidata/netcdf4-python.git)
* libhdf5-dev
* libnetcdf-dev
