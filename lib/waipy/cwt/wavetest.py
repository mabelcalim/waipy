#!/usr/bin/python
# -*- coding: latin-1 -*-
# WAVELET Torrence and Combo translate from Matlab to Python
# author: Mabel Calim Costa
# INPE
# 23/01/2013

"Baseado : Torrence e Combo"

# data from http://paos.colorado.edu/research/wavelets/software.html

import numpy as np
import pylab
from pylab import *
import matplotlib.pyplot as plt
import os
import sys
import netCDF4


def load_nc(file,var,dt,date1):
	"""OPEN ARCHIVE .NC
           file = archive.nc
	   var  = variable from archive.nc
	   dt 	= data sampling
	   date1= data intial time 
	"""
        f = netCDF4.Dataset(file,'r+')
        data = f.variables[var][:]
	n = len(data)
        time = np.arange(n)*dt + date1
        f.close()
        return data,time



def load_txt(archive,dt,date1):
	""" OPEN ARCHIVE .TXT/.DAT
	    archive = file.txt
	    dt 	    = data sampling
	    date1   = data initial time		
  	"""
        pkgs_pth = os.path.join(sys.prefix, 'lib', 'python'+sys.version[:3],'site-packages','wavelet','lib','wavelet','data','txt',archive)
        data = np.loadtxt(pkgs_pth)
        n = len(data)
        time = np.arange(n)*dt + date1  
	return data,time
#data,time = load_txt('sst_nino3.dat',0.25,1871)

	
def normalize(data):
        """NORMALIZE FUNCTION by - mean/sqrt(variance)
 	"""
	variance = np.var(data)
	data = (data - np.mean(data))/(np.sqrt(variance))
	return data
#data_norm = normalize(data)


def cwt (data, dt, pad, dj,s0,j1,lag1,param,mother): 
	"""CONTINUOUS WAVELET TRANSFORM
	   pad = 1         # pad the time series with zeroes (recommended)
	   dj = 0.25       # this will do 4 sub-octaves per octave
           s0 = 2*dt       # this says start at a scale of 6 months
           j1 = 7/dj       # this says do 7 powers-of-two with dj sub-octaves each
           lag1 = 0.72     # lag-1 autocorrelation for red noise background
           param = 6
           mother = 'Morlet'

        """
	variance = np.var(data)
	n = len(data)
	# Wavelet transform
	import lib_wavelet
	ondaleta,wave,period,scale,coi,f = lib_wavelet.wavelet(data,dt,mother,param)#,pad,dj,s0,j1,mother,param)
	power = np.abs(wave*wave)
	# Significance levels: (variance=1 for the normalized SST)
	signif,fft_theor = lib_wavelet.wave_signif(1.0,dt,scale,0,lag1,0.95,-1,mother,param)
	ones = np.ones((len(signif),n))		# expand signif --> ones (J+1)x(n)
	sig95 = [s * ones[1] for s in signif] 	# vstack signif concatenate
	sig95 = power/sig95			# where ratio > 1, power is significant
	# Global wavelet spectrum & significance levels:
	global_ws = variance*(np.sum(power.conj().transpose(),axis=0)/n)
	dof  = [n-s for s in scale]
	"""CAUTION - DEFAULT VALUES """
	global_signif,fft_theor = lib_wavelet.wave_signif(variance,dt,scale,1,lag1,0.95,dof,mother,param)
	# Daughter wavelet
	joint_wavelet = np.concatenate((np.fft.ifft(ondaleta)[np.ceil(n/2.):],np.fft.ifft(ondaleta)[np.ceil(n/2.):][::-1]), axis =1) #joint Wavelet
	imag_wavelet = np.concatenate((np.fft.ifft(ondaleta).imag[np.ceil(n/2.):],np.fft.ifft(ondaleta).imag[np.ceil(n/2.):][::-1]), axis =1)
	nw = np.size(joint_wavelet)  # daughter's number of points 
	# admissibility condition
	mean_wavelet = mean(joint_wavelet.real)
	mean_wavelet = np.ones(nw)*mean_wavelet
        result = {'ondaleta':ondaleta,'wave':wave ,'period': period,'scale':scale ,'coi':coi,'power': power,'sig95':sig95, 'global_ws':global_ws,
			'global_signif':global_signif,'joint_wavelet':joint_wavelet,'imag_wavelet':imag_wavelet,'nw':nw ,'mean_wavelet':mean_wavelet ,'dj':dj, 'j1':j1, 'dt':dt, 'fft':f }
	return result
 
#result = cwt(data_norm,0.25,1,0.25,2*0.25,7/0.25,0.72,6,'Morlet')

def fft(data):
        """FFT spectrum
        """
        n,X = len(data), np.fft.fft(data)
        sxx = (X * np.conj(X))/(n)
        f = -np.fft.fftfreq(n)[np.ceil(n/2.):]
        sxx = np.abs(sxx)
        sxx = sxx[np.ceil(n/2.):]
        return f,sxx

# ---------------------------
#           Ploting 
# ---------------------------
def levels(result,dtmin):
	""" Power levels
	    """
	#dtmin = result['power'].min()
	#dtmin = 0.0625
	dtmax = result['power'].max()
	lev = []
	for i in range(int(log2(dtmax/dtmin))):
        	dtmin =dtmin*2
        	lev.append(dtmin)
	return lev

def wavelet_plot(var,time,data,dtmin,result,impath):
	"""PLOT WAVELET TRANSFORM
	   var 	 = title name from data
	   time  = vector get in load function
	   data  = from normalize function
	   dtmin = minimum resolution :1 octave		 
	   result= dict from cwt function 
	"""
        from numpy import log2
        import numpy as np
	import wavetest
	fig = plt.figure(figsize=(15,10), dpi=100)
	ax = fig.add_subplot(411)
        # frequency limit
        lim = np.where(result['period'] == result['period'][-1]/2)[0][0]
        # ----------------------------------------------------------------------------------------------------------------#
	"""Plot time series """
	subplot(3,1,1)
	plot(time,data)
	xlabel('Time')
	ylabel('%s'%var)
	title('TIME SERIES')
	# ----------------------------------------------------------------------------------------------------------------#
	subplot(7,4,14)
	plot(range(-result['nw']/2,result['nw']/2),result['joint_wavelet'],'k')
	plot(range(-result['nw']/2,result['nw']/2),result['imag_wavelet'],'--k')
	plot(range(-result['nw']/2,result['nw']/2),result['mean_wavelet'],'g')
	xlim(-50,50)
	ylim(-0.2,0.2)
	text(20,0.10,'Morlet')
	title('$\psi$ (t/s)')
	# ----------------------------------------------------------------------------------------------------------------#
	subplot(7,4,15)
	plot(result['ondaleta'],'k') 
	title('$\psi^-$  ')
	# ----------------------------------------------------------------------------------------------------------------#
	""" Contour plot wavelet power spectrum """
	subplot(3,2,5)
	lev = wavetest.levels(result,dtmin)									# call levels - power levels
	contourf(time, np.log2(result['period']),np.log2(result['power']),np.log2(lev))
	colorbar()
	ax = gca()											# handle to the current axes
	ax.set_ylim(ax.get_ylim()[::-1])								# reverse plot along y axis
	yt =  range(int(np.log2(result['period'][0])),int(np.log2(result['period'][-1])+1)) 		# create the vector of periods
	Yticks = [float(math.pow(2,p)) for p in yt]							# make 2^periods
	yticks(yt, map(str,Yticks))
	xlim(time[0],time[-1])                  							# date on the x axis 
	# 95% significance contour, levels at -99 (fake) and 1 (95% signif)
	contour(time,np.log2(result['period']),result['sig95'], [-99,1],linewidths= 2)
	# cone-of-influence , anything "below"is dubious
	plot(time,np.log2(result['coi']),'k')
	#ax.fill(time,np.log2(result['coi']),'k', alpha =0.3, hatch = '/') 
	ax.fill_between(time,np.log2(result['coi']),int(np.log2(result['period'][-1])+1), alpha =0.5, hatch = '/')
	xlabel('Time')
	ylabel('Period')
	title('Wavelet Power Spectrum')
	# ----------------------------------------------------------------------------------------------------------------#
	""" Plot global wavelet spectrum"""
	subplot(3,2,6)
        f,sxx = wavetest.fft(data)
        #yt_fft = np.linspace(result['period'][0],result['period'][lim],len(f))
        plot(sxx,np.log2(1/f*result['dt']),'gray', label = 'Fourier spectrum')
        yticks(yt, map(str,Yticks))
	plot(result['global_ws'], np.log2(result['period']),'b', label = 'Wavelet spectrum')
	plot(result['global_signif'],np.log2(result['period']),'r--')
	ax = gca()                              # handle to the current axes
	ax.set_ylim(ax.get_ylim()[::-1])        # reverse plot along y axis
	yticks(yt, map(str,Yticks))
	##xlim(0,1.25*np.max(global_ws))
	##xlabel('Power (units of %var $^{2})$'%var)
	xlabel('Power')
	title('Global Wavelet Spectrum')
        #legend = ax.legend(loc='lower right', shadow=True)
	#plt.savefig('%var.png'%var)
        plt.savefig(os.path.join(impath,'%s.png'%var))
	pylab.show()
        return
#wavelet_plot('SST_NINO3',time,data,0.03125,result,impath)
