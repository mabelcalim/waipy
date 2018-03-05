waipy
=====
This guide includes a Continuous Wavelet Transform (CWT), significance  tests
from based on Torrence and Compo (1998)
![ScreenShot](https://github.com/mabelcalim/waipy/blob/master/Sine.png)


and Cross Wavelet Analysis  (CWA) based on Maraun and Kurths(2004).

Cross Power Wavelet Analysis:
    ![ScreenShot](https://github.com/mabelcalim/waipy/blob/master/Waipy%20Examples%20/Jupyter/figs/CrossWavelet_noise_example.png)

Wavelet Coherence Analysis:
    ![ScreenShot](https://github.com/mabelcalim/waipy/blob/master/Waipy%20Examples%20/Jupyter/figs/CohereWavelet_noise_example.png)


Installation
============

First steps:

    cd  PATH_TO_SAVE_WAIPY_IN_YOUR_HOME   
    git clone https://github.com/mabelcalim/waipy.git
    cd waipy/

A local install can be done using the provided setup.py file:

    python3 setup.py install

The installation path can be changed using the **--prefix** switch, e.g.:

    python setup.py install --prefix $HOME/inst

Make sure to add the corresponding paths to your ``$PATH`` and ``$PYTHONPATH``
environment variables. Alternatively, if the **--user** switch can be used:

    export PYTHONUSERBASE=$HOME/inst/pip_installs
    export PYTHONPATH=$HOME/inst/pip_installs/lib/python2.7/\
        site-packages/:$PYTHONPATH
    export PATH=$HOME/inst/pip_installs/bin:$PATH
    python setup.py install --user

If you plan on modifying the code, use the **develop** target in combination
with the **--user** swich:

    export PYTHONUSERBASE=$HOME/inst/pip_installs
    export PYTHONPATH=$HOME/inst/pip_installs/lib/python2.7/\
        site-packages/:$PYTHONPATH
    export PATH=$HOME/inst/pip_installs/bin:$PATH
    python setup.py develop --user

The first three lines should also be included in the **$HOME/.bashrc** file.

Requirements
============

The following packages are required (tested on Debian Wheezy):

* python-matplotlib
* python-numpy
* libnetcdf-dev
* libhdf5-dev
* netcdf4 (https://github.com/Unidata/netcdf4-python.git)


Installing requirements for python3
============

* pip3 install numpy
* pip3 install matplotlib
* pip3 install netCDF4
* pip3 install pandas
* pip3 install scipy

As a tip, try to use jupyter !
* pip3 install jupyter


Examples
============

* https://github.com/mabelcalim/waipy/tree/master/Waipy%20Examples%20/Jupyter

Use waipy without install it!
============
Check the prêt-à-porter.ipynb examples!!!
* https://github.com/mabelcalim/waipy/blob/master/Waipy%20Examples%20/waipy_pr%C3%AAt-%C3%A0-porter.ipynb
* https://github.com/mabelcalim/waipy/blob/master/Waipy%20Examples%20/cwa_pr%C3%AAt-%C3%A0-porter.ipynb


Support Group
==============

https://groups.google.com/forum/?hl=en#!forum/waipy-users-support
