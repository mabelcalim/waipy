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


# load data -- transf .dat --> array
sst = np.loadtxt('can_filt.txt') 
#sst = np.loadtxt('sst_nino3.dat')

""" Normalize by standard deviation """
variance = np.var(sst)
sst = (sst - np.mean(sst))/np.sqrt(variance)

n = len(sst)
dt = 0.25
time = np.arange(n)*dt + 1871
# plt.xlim(1870,200)
pad = 1  	# pad the time series with zeroes (recommended)
dj = 0.25	# this will do 4 sub-octaves per octave
s0 = 2*dt	# this says start at a scale of 6 months
j1 = 7/dj	# this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72	# lag-1 autocorrelation for red noise background
param = 6
mother = 'Morlet'
#param = 2
#mother = 'DOG'

# Wavelet transform
import lib_wavelet
ondaleta,wave,period,scale,coi = lib_wavelet.wavelet(sst,dt,mother,param)#,pad,dj,s0,j1,mother,param)
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




# ---------------------------
#           Ploting 
# ---------------------------
fig = plt.figure(figsize=(15,10), dpi=100)
ax = fig.add_subplot(411)
"""Plot time series """
subplot(3,1,1)
plot(time,sst)
xlabel('Time (year)')
ylabel('NINO3 SST (degC)')
title('a) NINO3 Sea Surface Temperature (seasonal)')

subplot(7,4,14)
plot(range(-nw/2,nw/2),joint_wavelet,'k')
plot(range(-nw/2,nw/2),imag_wavelet,'--k')
plot(range(-nw/2,nw/2),mean_wavelet,'g')
xlim(-50,50)
ylim(-0.2,0.2)
text(20,0.10,'Morlet')
title('$\psi$ (t/s)')

subplot(7,4,15)
plot(ondaleta,'k') 
title('$\psi^-$  ')

""" Contour plot wavelet power spectrum """
subplot(3,2,5)
levels = [0.0625,0.125,0.25,0.5,1,2,4,8,16]
contour(time, np.log2(period),np.log2(power),np.log2(levels))
ax = gca()							# handle to the current axes
ax.set_ylim(ax.get_ylim()[::-1])				# reverse plot along y axis
yt =  range(int(np.log2(period[0])),int(np.log2(period[-1]))+1) # create the vector of periods
Yticks = [int(math.pow(2,p)) for p in yt]			# make 2^periods
yticks(yt, map(str,Yticks))
xlim(time[0],time[-1])                  			# date on the x axis 
# 95% significance contour, levels at -99 (fake) and 1 (95% signif)
contour(time,np.log2(period),sig95, [-99,1],linewidths= 2)
# cone-of-influence , anything "below"is dubious
#coi = np.ma.masked_where(coi <= np.log2(power),coi)
plot(time,np.log2(coi),'k')


""" CAUTION  - shaded coi area"""
# shaded coi area
#ax.fill_between(time,np.log2(period),np.log2(coi), facecolor='gray', interpolate=True)

xlabel('Time (year)')
ylabel('Period (years)')
title('b) NINO3 SST Wavelet Power Spectrum')


""" Plot global wavelet spectrum"""
subplot(3,2,6)
plot(global_ws, np.log2(period))
plot(global_signif,np.log2(period),'r--')
ax = gca()                              # handle to the current axes
ax.set_ylim(ax.get_ylim()[::-1])        # reverse plot along y axis
yticks(yt, map(str,Yticks))
#xlim(0,1.25*np.max(global_ws))
xlabel('Power (deg C$^{2})$')
title('c) Global Wavelet Spectrum')

#plt.savefig('wavelet_nino.png')


pylab.show()
