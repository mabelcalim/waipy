# !/usr/bin/python
# -*- coding: latin-1 -*-
# WAVELET Torrence and Combo translate from Matlab to Python
# author: Mabel Calim Costa
# INPE
# 23/01/2013

#Baseado : Torrence e Combo

# data from http://paos.colorado.edu/research/wavelets/software.html

import numpy as np
import pylab
from pylab import *
import matplotlib.pyplot as plt
import os
import sys
import netCDF4


def load_nc(file, var, dt, date1):
    """
    Opens and reads an archive .txt with data only (without dates).
    """
    f = netCDF4.Dataset(file, 'r+')
    data = f.variables[var][:]
    n = len(data)
    time = np.arange(n) * dt + date1
    f.close()
    return data, time


def load_txt(archive, dt, date1):
    """
    Opens and reads an archive .txt with data only (without dates).

    :Parameters:
      |  **file:** filename
      |  **dt:** A number. Time-step of the vector. Example: Hourly, daily, monthly, etc.
      |  **date1:** A number. The initial time of the data. Example: 1985.
    :Returns:
      |  **data:** array_like
      |  **time:** array_like

    :Notes:
      |  This function is linked to data/txt directory, so, if you have any file extension .txt put it in the following folder:

    ::

        /lib/wavelet/data/txt

    Example
    ::

        >> dt = 0.25
        >> date1 = 1871
        # Test data = sst_nino3.dat is already in the package!
        >> data, time = load_txt('sst_nino3.dat', dt, date1)
    """

    filename = os.path.join(sys.prefix, 'lib', 'python' + sys.version[
                            :3], 'site-packages', 'wavelet', 'lib',
                            'wavelet', 'data', 'txt', archive)
    if(not os.path.isfile(filename) and os.path.isfile(archive)):
        filename = archive
    else:
        raise IOError(
            'File {0} not found either here {1} or here {1}'.format(filename,
                                                                    archive))
    data = np.loadtxt(filename)
    n = len(data)
    time = np.arange(n) * dt + date1
    return data, time


def normalize(data):
    """
    Normalize by Standard Score - mean = 0 ; Variance = 1.
    
    :Parameters:
      |  **Data, the loaded data.**
    :Returns:
      |  **normalized_data**

    :Notes:
      |  You can skip this function if it the normalization is not necessary (e.g. EOF data).
    """
    variance = np.var(data)
    data = (data - np.mean(data)) / (np.sqrt(variance))
    return data


def cwt(data, dt=1, pad=1, dj=0.25, s0=2, j1=7/0.25, lag1=0.72, param=6, mother='Morlet', name='Name of Time series'):
    """
    Continuous wavelet transform from data. Wavelet params can be modified as you wish.

    :Parameters:
      |  **Data:** array_like. Raw of data or normalized data.
      |  **pad:** 0 or 1 Pad the time series with zeroes to next pow of two length (recommended).
      |  **dt:** Time-step of the vector. Example: Hourly, daily, monthly, etc...
      |  **dj:** Divide octave in sub-octaves. If dj = 0.25 this will do 4 sub-octaves per octave.
      |  **s0:** The maximum frequency resolution. If it is an annual data, s0 = 2*dt say start at a scale of 6 months.
      |  **j1:** Divide the power-of-teo with dj sub-octaves each.
      |  **lag1:** Lag-1 autocoorelation for red noise background.
      |  **param:** The mother wavelet nondimensional time-parameter, depends on wavelet, for other wavelet than Morlet, see Torrence and Compo 1998.
      |  **mother:** string. The mother wavelet funtion. Can be 'Morlet', 'PAUL', 'DOG'.
      |  **name:** Name of Time Series or your Plot.

    :Returns:
      |  **result:** as dict. Returns all parameters for plot

    Example
    ::

        >> dt = 0.25
        >> date1 = 1871
        # Test data = sst_nino3.dat is already in the package!
        >> data,n,time = load_txt('sst_nino3.dat',dt,date1)
        # This normalize by variance
        >> data_norm, variance = normalize(data)
        # Continuous wavelet transform
        >> result = cwt(data_norm,0.25,variance,n,1,0.25,2*0.25,7/0.25,0.72,6,'Morlet')
    """

    import lib_wavelet

    variance = np.var(data)
    n = len(data)
    # Wavelet transform
    ondaleta, wave, period, scale, coi, f = lib_wavelet.wavelet(
        data, dt, param, dj, s0, j1, mother)
    # wave = np.array(wave)
    power = (np.abs(wave) ** 2)
    # Significance levels: (variance=1 for the normalized SST)
    signif, fft_theor = lib_wavelet.wave_signif(
        1.0, dt, scale, 0, lag1, 0.95, -1, mother, param)
    ones = np.ones((len(signif), n))  # expand signif --> ones (J+1)x(n)
    sig95 = [s * ones[1] for s in signif]   # vstack signif concatenate
    sig95 = power / sig95  # where ratio > 1, power is significant
    # Global wavelet spectrum & significance levels:
    global_ws = variance * (np.sum(power.conj().transpose(), axis=0) / n)
    dof = [n - s for s in scale]
    """CAUTION - DEFAULT VALUES """
    global_signif, fft_theor = lib_wavelet.wave_signif(
        variance, dt, scale, 1, lag1, 0.95, dof, mother, param)
    # Daughter wavelet
    joint_wavelet = np.concatenate((np.fft.ifft(ondaleta)[np.ceil(
        n / 2.):], np.fft.ifft(ondaleta)[np.ceil(n / 2.):][::-1]), axis=1)
    imag_wavelet = np.concatenate((np.fft.ifft(ondaleta).imag[np.ceil(
        n / 2.):], np.fft.ifft(ondaleta).imag[np.ceil(n / 2.):][::-1]), axis=1)
    nw = np.size(joint_wavelet)  # daughter's number of points
    # admissibility condition
    mean_wavelet = mean(joint_wavelet.real)
    mean_wavelet = np.ones(nw) * mean_wavelet
    result = {'ondaleta': ondaleta, 'wave': wave, 'period': period,
              'scale': scale, 'coi': coi, 'power': power, 'sig95': sig95,
              'global_ws': global_ws, 'global_signif': global_signif,
              'joint_wavelet': joint_wavelet, 'imag_wavelet': imag_wavelet,
              'nw': nw, 'mean_wavelet': mean_wavelet, 'dj': dj, 'j1': j1,
              'dt': dt, 'fft': f, 'mother': mother, 'data': data, 'name': name}
    return result

# result = cwt(data_norm,0.25,1,0.25,2*0.25,7/0.25,0.72,6,'Morlet')


def fft(data):
    """FFT spectrum
    """
    n = len(data)
    X = np.fft.fft(data)
    sxx = ((X * np.conj(X)) / (n))
    f = -np.fft.fftfreq(n)[np.ceil(n / 2.):]
    sxx = np.abs(sxx)
    sxx = sxx[np.ceil(n / 2.):]
    return f, sxx

# ---------------------------
#           Ploting
# ---------------------------


def levels(result, dtmin):
    """
    Power levels
    """

    dtmax = result['power'].max()
    lev = []
    for i in range(int(log2(dtmax / dtmin))):
        dtmin = dtmin * 2
        lev.append(dtmin)
    return lev


def wavelet_plot(var, time, data, dtmin, result):
    """
    PLOT WAVELET TRANSFORM
    var = title name from data
    time  = vector get in load function
    data  = from normalize function
    dtmin = minimum resolution :1 octave
    result = dict from cwt function
    """

    from numpy import log2
    import numpy as np
    import wavetest
    import matplotlib
    import matplotlib.gridspec as gridspec

    # frequency limit
    # print result['period']
    # lim = np.where(result['period'] == result['period'][-1]/2)[0][0]
    """Plot time series """

    fig = plt.figure(figsize=(15, 10), dpi=300)

    gs1 = gridspec.GridSpec(4, 3)
    gs1.update(left=0.05, right=0.7, wspace=0.5, hspace=0)
    ax1 = plt.subplot(gs1[0, :])
    ax1 = pylab.gca()
    ax1.xaxis.set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs1[1:4, :])#, axisbg='#C0C0C0')

    gs2 = gridspec.GridSpec(4, 1)
    gs2.update(left=0.7, right=0.98, hspace=0)
    ax5 = plt.subplot(gs2[1:4, 0], sharey=ax2)
    plt.setp(ax5.get_yticklabels(), visible=False)

    gs3 = gridspec.GridSpec(6, 1)
    gs3.update(left=0.77, top=0.86, right=0.98, hspace=0.6, wspace=0.01)
    ax3 = plt.subplot(gs3[0, 0])

    # ----------------------------------------------------------------------------------------------------------------#
    ax1.plot(time, data)
    ax1.axis('tight')
    ax1.set_ylabel('SP [mV]', fontsize=15)
    ax1.set_title('%s' % var, fontsize=17)
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax1.grid(True)
    ax1.xaxis.set_visible(False)
    # ----------------------------------------------------------------------------------------------------------------#
    ax3.plot(range(-result['nw'] / 2, result['nw'] / 2),
             result['joint_wavelet'], 'k', label='Real part')
    ax3.plot(range(-result['nw'] / 2, result['nw'] / 2),
             result['imag_wavelet'], '--k', label='Imag part')
    ax3.plot(range(-result['nw'] / 2, result['nw'] / 2),
             result['mean_wavelet'], 'g', label='Mean')
    # ax3.axis('tight')
    ax3.set_xlim(-40, 40)
    # ax3.set_ylim(-0.3,0.3)
    # ax3.set_ylim([np.min(result['joint_wavelet']),np.max(result['joint_wavelet'])])
    ax3.set_xlabel('Time', fontsize=10)
    ax3.set_ylabel('Amplitude', fontsize=10)
    ax3.set_title('$\psi$ (t/s) {0} in time domain'.format(result['mother']))
    # ----------------------------------------------------------------------------------------------------------------#
    # ax4.plot(result['ondaleta'],'k')
    # ax4.set_xlabel('Frequency', fontsize=10)
    # ax4.set_ylabel('Amplitude', fontsize=10)
    # ax4.set_title('$\psi^-$  Frequency domain', fontsize=13)
    # ----------------------------------------------------------------------------------------------------------------#
    """ Contour plot wavelet power spectrum """
    lev = wavetest.levels(result, dtmin)
    pc = ax2.contourf(
        time, np.log2(result['period']),
        np.log2(result['power']), np.log2(lev))
    # 95% significance contour, levels at -99 (fake) and 1 (95% signif)
    pc2 = ax2.contour(
        time, np.log2(result['period']), result['sig95'],
        [-99, 1], linewidths=2)
    ax2.plot(time, np.log2(result['coi']), 'k')
    # cone-of-influence , anything "below"is dubious
    ax2.fill_between(time, np.log2(result['coi']), int(
        np.log2(result['period'][-1]) + 1), alpha=0.5, hatch='/')
    position = fig.add_axes([0.07, 0.07, 0.6, 0.01])
    cbar = plt.colorbar(pc, cax=position, orientation='horizontal')
    cbar.set_label('Power')
    yt = range(int(np.log2(result['period'][0])), int(
        np.log2(result['period'][-1]) + 1))  # create the vector of periods
    Yticks = [float(math.pow(2, p)) for p in yt]  # make 2^periods
    # Yticks = [int(i) for i in Yticks]
    ax2.set_yticks(yt)
    ax2.set_yticklabels(Yticks)
    ax2.set_ylim(ymin=(np.log2(np.min(result['period']))), ymax=(
        np.log2(np.max(result['period']))))
    ax2.set_ylim(ax2.get_ylim()[::-1])
    ax2.set_xlabel('Time', fontsize=12)
    ax2.set_ylabel('Period', fontsize=12)
    ax2.axhline(y=10.5, xmin=0, xmax=1, linewidth=2, color='k')
    ax2.axhline(y=13.3, xmin=0, xmax=1, linewidth=2, color='k')
    # ----------------------------------------------------------------------------------------------------------------#
    """ Plot global wavelet spectrum """
    f, sxx = wavetest.fft(data)
    ax5.plot(
        sxx, np.log2(1 / f * result['dt']), 'gray', label='Fourier spectrum')
    ax5.plot(result['global_ws'], np.log2(
        result['period']), 'b', label='Wavelet spectrum')
    ax5.plot(result['global_signif'], np.log2(
        result['period']), 'r--', label='95% confidence spectrum')
    ax5.legend(loc=0)
    ax5.set_xlim(0, 1.25 * np.max(result['global_ws']))
    ax5.set_xlabel('Power', fontsize=10)
    ax5.set_title('Global Wavelet Spectrum', fontsize=12)
    # save fig
    plt.savefig('%s.png' % var, dpi=300)
