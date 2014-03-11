# -*- coding: utf-8 -*-
#!/usr/bin/python
# Cross Wavelet Analysis (CWA) based on Maraun and Kurths(2004).
# http://www.nonlin-processes-geophys.net/11/505/2004/npg-11-505-2004.pdf
# author: Mabel Calim Costa
# INPE
# 23/01/2013

"""
Created on Mon Jun 17 2013

@author: Mabel Calim Costa
"""
import numpy as np
import pylab
from pylab import *
import matplotlib.pyplot as plt



def cross_wavelet (wave1,wave2):
    	""" Computes the cross wavelet analysis.
		wave1 = result['wave'] time serie 1
        	wave2 = result['wave'] time serie 2
        A normalized time and scale resolved measure for the relationship 
        between two time series x1(ti) and x2(ti) is the wavelet coherency (WCO),
        which is defined as the amplitude of the WCS(wavelet cross spectrum) 
        normalized to the two single WPS(wavelet power spectrum) (Maraun and Kurths,2004).
         	WCOi(s) = |WCSi(s) |/ (WPSi1(s) WPSi2(s)) ˆ1/2 
        _____________________________________________________________________
        Inputs:  
		wave1 - wavelet transform of time series x1
                wave2 - wavelet transform of time series x2
        Outputs:
		cohere - wavelet coherency (WCO)
 	Call function:
		cohere = cross_wavelet(wave1,wave2)
    	""" 
    	cross_power = np.abs(wave1.real*wave2.imag)
   	coherence = np.sqrt(cross_power*cross_power)/np.sqrt(np.abs(wave1.real*wave1.imag)* np.abs(wave2.real*wave2.imag))
   	return cross_power, coherence


def plot_cross (cross_power,time,result):
	""" PLOT CROSS POWER
	    cross_power = from cross_wavelet function
	    coherence   = from cross_wavelet function
	"""
	fig = plt.figure(figsize=(15,10), dpi=100)
	CS = plt.contourf(time,np.log2(result['period']),cross_power)
	cbar = plt.colorbar(CS)
	ax = gca()                                                      # handle to the current axes
	ax.set_ylim(ax.get_ylim()[::-1])                                # reverse plot along y axis
	yt =  range(int(np.log2(result['period'][0])),int(np.log2(result['period'][-1]))+1) # create the vector of periods
	xt = range(int(time[0]),int(time[-1]))
	Yticks = [(math.pow(2,p)) for p in yt]         		        # make 2^periods
	yticks(yt, map(str,Yticks))
	xticks(xt, map(str,xt))
	xlim(time[0],time[-1])                                          # date on the x axis 
	xlabel('Time')
	ylabel('Period')
	title('Cross Power')
        plot(time,np.log2(result['coi']),'k')
        ax.fill_between(time,np.log2(result['coi']),int(np.log2(result['period'][-1])*2), alpha =0.5, hatch = '/')
	pylab.show()
	return

def plot_cohere (coherence,time,result):
	"""PLOT COHERENCE
	   coherence   =  from cross_wavelet function
	   time	       =  time vector from load function
	   result      =  dict from cwt function
        """
        fig = plt.figure(figsize=(15,10), dpi=100)
        lev = list(np.linspace(0,1.0, 6))
        CS = plt.contourf(time,np.log2(result['period']),coherence,lev)
	cbar = plt.colorbar(CS)
        ax = gca()                                                      # handle to the current axes
        ax.set_ylim(ax.get_ylim()[::-1])                                # reverse plot along y axis
        yt =  range(int(np.log2(result['period'][0])),int(np.log2(result['period'][-1]))+1) # create the vector of periods
	xt = range(int(time[0]),int(time[-1]))
        Yticks = [(math.pow(2,p)) for p in yt]                          # make 2^periods
        yticks(yt, map(str,Yticks))
	xticks(xt, map(str,xt))
        xlim(time[0],time[-1])                                          # date on the x axis 
        xlabel('Time')
        ylabel('Period')
        title('Coherence')
        plot(time,np.log2(result['coi']),'k')
        ax.fill_between(time,np.log2(result['coi']),int(np.log2(result['period'][-1])*2), alpha =0.5, hatch = '/')
        pylab.show()
        return

