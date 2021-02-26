# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:03 2013

@author: Mabel Calim Costa
"""

import os
from setuptools import setup
#from distutils.core import setup


for line in open('lib/waipy/__init__.py').readlines():
    if line.startswith('__version__'):
        exec(line.strip())

setup(
    name = "waipy",
    description = ("Wavelet Analysis in Python"),
    version=__version__,
    author='Mabel Calim Costa',
    author_email='mabelcalim@gmail.com',
    #url='https://wavelet-analysis.readthedocs.org/en/latest/index.html',
    url = 'https://bitbucket.org/mabel/waipy/overview',
    long_description="""
This guide includes a Continuous Wavelet Transform (CWT), significance
tests from based on Torrence and Compo (1998) and Cross Wavelet Analysis
(CWA) based on Maraun and Kurths(2004).""",
    packages=['waipy', 'waipy.cwt', 'waipy.cwa' ],
    package_dir={'':'lib'},
    classifiers=['License :: OSI Approved :: BSD License'],
    #install_requires=['numpy', 'scipy','pandas', 'matplotlib'],
    extras_require= {
        'all': ["netCDF4", "jupyter"],
        'load_netcdf': ["netCDF4"],
        'jupyer': ["jupyter"],
    },
)
