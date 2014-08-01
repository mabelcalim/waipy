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
import cmath
import pandas as pd

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

    cross_power = wave1*wave2.conjugate()
    cross_power1 = np.abs(wave1*wave2.conjugate())
    WPS1 = np.abs(wave1*wave1.conjugate())
    WPS2 = np.abs(wave2*wave2.conjugate())
    #   	coherence = np.sqrt(cross_power*cross_power)/np.sqrt(np.abs(wave1.real*wave1.imag)* np.abs(wave2.real*wave2.imag))
    coherence = cross_power/(np.sqrt(WPS1*WPS2))
    phase_angle = np.angle(cross_power)#,deg=True)
    return cross_power1, coherence, phase_angle


def plot_cross (var,cross_power,phase_angle, time, result, result1):
    """ PLOT CROSS POWER
    cross_power = from cross_wavelet function
    coherence   = from cross_wavelet function
    """
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(15,10),dpi=300) 
    #set plot grid
    gs1 = gridspec.GridSpec(4, 3)
    gs1.update(left=0.05)
    ax1 = plt.subplot(gs1[0, :])
    ax1 = pylab.gca()
    ax1.xaxis.set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs1[1:4, :], axisbg='#C0C0C0')
    #plot timeseries
    ax1.plot(time,result['data'])
    ax1.set_ylabel('SP [mV]', fontsize=15)
    ax3 = ax1.twinx()
    ax3.plot(time,result1['data'], color='c')
    ax3.set_ylabel(result1['name'])
    ax1.set_title('%s'%var, fontsize=17)
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax1.grid(True)
    ax1.xaxis.set_visible(False)

    phs_dt = round(len(time)/20)
    tidx = np.arange(np.max(np.floor(phs_dt/2)),len(time),phs_dt)
    tidx = [int(i) for i in tidx]
    tidx = np.array(tidx)
    #print tidx
    phs_dp = round(len(result['period'])/20) 
    pidx = np.arange(np.max(np.floor(phs_dp/2)),len(result['period']),phs_dp)
    pidx = [int(i) for i in pidx]
    pidx = np.array(pidx)
    X,Y = meshgrid(time.astype(np.int64)[tidx],np.log2(result['period'][pidx]))
    phase_angle1 = phase_angle[:,tidx]
    phase_angle1 = phase_angle1[pidx,:]
    cA = np.exp(1j * phase_angle1)
    U = np.real(cA)
    V = np.imag(cA)
    #X = pd.DatetimeIndex(X)
    #print V.shape
    #print U.shape
    #print 
    #print Y.shape
    #plot cross
    ax4 = ax2.twiny()
    ax4.xaxis.set_visible(False)
    #ax4.set_xlim(0.9,4.4)
    CS = ax2.contourf(time, np.log2(result['period']), cross_power)
    ax2.plot(time, np.log2(result['coi']), 'k')# cone-of-influence , anything "below"is dubious
    ax2.fill_between(time,np.log2(result['coi']), int(np.log2(result['period'][-1])+1), alpha =0.5, hatch = '/')
    position=fig.add_axes([0.07,0.07,0.6,0.01])  ##
    cbar = plt.colorbar(CS, cax=position,orientation='horizontal')#,norm=normal)#,shrink=0.5,pad=0.08)
    cbar.set_label('Power')
    Q = ax4.quiver(X.astype(np.int64), Y, U, V)#,angles='uv',linewidth=1)
    ax4.axis('tight')
    yt =  range(int(np.log2(result['period'][0])), int(np.log2(result['period'][-1])+1)) # create the vector of periods
    Yticks = [float(math.pow(2,p)) for p in yt] # make 2^periods
    Yticks = [int(i) for i in Yticks] 
    ax2.set_yticks(yt)
    ax2.set_yticklabels(Yticks)
    ax2.set_ylim(ymin=(np.log2(result['period'][0])), ymax=(np.log2(result['period'][-1])))
    ax2.invert_yaxis()
    ax2.set_xlabel('Time', fontsize=15)
    ax2.set_ylabel('Period (Minutes)', fontsize=15)
    ax2.axhline(y=10.5, xmin=0, xmax=1, linewidth=2, color='k')
    ax2.axhline(y=13.3, xmin=0, xmax=1, linewidth=2, color='k')
    ax2.set_title('Cross Power')






    plt.savefig('Cross Power {0} vs {1}'.format(result['name'],result1['name']),dpi=300)

def plot_cohere (var,coherence, time, result, result1):
    """
    PLOT COHERENCE
    coherence   =  from cross_wavelet function
    time       =  time vector from load function
    result      =  dict from cwt function
    """
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(15,14),dpi=300) 
    #set plot grid
    gs1 = gridspec.GridSpec(4, 3)
    gs1.update(left=0.05)
    ax1 = plt.subplot(gs1[0, :])
    ax1 = pylab.gca()
    ax1.xaxis.set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs1[1:4, :], axisbg='#C0C0C0')

    #plot timeseries
    ax1.plot(time,result['data'])
    ax1.set_ylabel('SP [mV]', fontsize=15)
    ax3 = ax1.twinx()
    ax3.plot(time,result1['data'], color='c')
    ax3.set_ylabel(result1['name'])
    ax1.set_title('%s'%var, fontsize=17)
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax1.grid(True)
    ax1.xaxis.set_visible(False)

    #fig = plt.figure(figsize=(15,10), dpi=100)
    lev = list(np.linspace(0,1.0, 6))
    CS = ax2.contourf(time,np.log2(result['period']),coherence,lev)
    ax2.plot(time,np.log2(result['coi']),'k')
    ax2.fill_between(time,np.log2(result['coi']),int(np.log2(result['period'][-1])*2), alpha =0.5, hatch = '/')
    position=fig.add_axes([0.07,0.07,0.6,0.01])  ##
    cbar = plt.colorbar(CS, cax=position,orientation='horizontal')#,norm=normal)#,shrink=0.5,pad=0.08)
    #cbar.set_label('Power')
    yt =  range(int(np.log2(result['period'][0])), int(np.log2(result['period'][-1])+1)) # create the vector of periods
    Yticks = [float(math.pow(2,p)) for p in yt] # make 2^periods
    Yticks = [int(i) for i in Yticks] 
    ax2.set_yticks(yt)
    ax2.set_yticklabels(Yticks)
    ax2.set_ylim(ymin=(np.log2(result['period'][0])), ymax=(np.log2(result['period'][-1])))
    ax2.invert_yaxis()
    ax2.set_xlabel('Time', fontsize=15)
    ax2.set_ylabel('Period (Minutes)', fontsize=15)
    #ax2.axhline(y=10.5, xmin=0, xmax=1, linewidth=2, color='k')
    #ax2.axhline(y=13.3, xmin=0, xmax=1, linewidth=2, color='k')
    ax2.set_title('Coherence')
    plt.savefig('Coherence {0} vs {1}'.format(result['name'],result1['name']),dpi=300)
    return  
#cbar = plt.colorbar(CS)
    #ax = gca()                                                      # handle to the current axes
    #ax.set_ylim(ax.get_ylim()[::-1])                                # reverse plot along y axis
    #yt =  range(int(np.log2(result['period'][0])),int(np.log2(result['period'][-1]))+1) # create the vector of periods
    #xt = range(int(time[0]),int(time[-1]))
    #Yticks = [(math.pow(2,p)) for p in yt]                          # make 2^periods
    #yticks(yt, map(str,Yticks))
    #xticks(xt, map(str,xt),rotation=45)
    #xlim(time[0],time[-1])                                          # date on the x axis 
    #xlabel('Time')
    #ylabel('Period')
    #title('Coherence')
    


