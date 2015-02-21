#!/usr/bin/env Python

# This module reads a configuration file containing 
# the wfdisc tuples and waveform files in CSS 
# four byte integer s4 big endian format and:
#
# 1 - Selects a portion of the waveform as master waveform for crosscorrelation
# 2 - Displays the master waveform 
# 2 - Displays the waveforms, crosscorrelation of the three instruments within a 
#     triad
# 3 - Displays the picks on the crosscorrelation when they are larger than a threshold
#    
# 
# Developped February 2013 - Ronan Le Bras - Gydatos LLC
#

# import all necessary modules

import struct
import math
from   array import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 
from   matplotlib.widgets import Button, Slider
import ConfigParser
import string
import time
import calendar
import scipy.signal
import numpy as np
import scipy.io.matlab.byteordercodes
import pylab
#import pyaudio
import wave
import glob
import os
import subprocess
import gydatos as gy

epsilon = 0.0001



# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('Waveforms.cfg')

verbose = config.getint('control_parameters', 'verbose')

#wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

envelope = config.getboolean('Filter','envelope')
order = config.getint('Filter','order')
low = config.get('Filter','low')
high = config.get('Filter','high')

start_time = config.get('time_interval','start_time')
end_time = config.get('time_interval','end_time')
time_format = config.get('time_interval','format') 
time_interval = config.getfloat('time_interval','time_interval') 
harmonic_freq = config.getfloat('time_interval','harmonic_frequency') 
harmonic_length = config.getfloat('time_interval','harmonic_length') 


# station = config.get('time_interval','Station')

triad = config.get('Correlation','Triad')

station1 = triad+'1'
station2 = triad+'2'
station3 = triad+'3'

wf1=array('f',[])
wf2=array('f',[])
wf3=array('f',[])
wfm=array('f',[])

#sta_master   = config.get('Correlation','station_master')
#stime_master = config.get('Correlation','start_master')
#etime_master = config.get('Correlation','end_master')
#halfwidth = config.get('Azimuth','halfwidth')
srate=gy.ReadSrate(wfdisc, station1, verbose)
print ("Sampling rate:", srate,'\n')
if verbose > 2:
	print 'Order =', order, 'Low =', low, 'High =', high


# Parse time variables and get start and end epoch time

stime = time.strptime(start_time, time_format)
etime = time.strptime(end_time, time_format)

start_e = calendar.timegm(stime)
end_e   = calendar.timegm(etime)

for i in range(1+int(time_interval*srate)):
		wf1.append(0.)
for i in range(1+int(time_interval*srate)):
		wf2.append(0.)
for i in range(1+int(time_interval*srate)):
		wf3.append(0.)	
	
if verbose > 1:
	print "Start time:", stime,'\n'
	print "Epoch start:", start_e
	print "End time:", etime,'\n'
	print "Epoch end:", end_e

samprate, start_wf1, num_samp1 = gy.ReadWfm(wf1, wfdisc, start_e, end_e, station1, verbose)
samprate, start_wf2, num_samp2 = gy.ReadWfm(wf2, wfdisc, start_e, end_e, station2, verbose)
samprate, start_wf3, num_samp3 = gy.ReadWfm(wf3, wfdisc, start_e, end_e, station3, verbose)

#samprate, null, wfm = gy.ReadWfm(wfdisc, start_master_e, end_master_e, sta_master, verbose)


# Design digital IIR filter

Wn = array('f',[])
Wn.append(float(low))
Wn.append(float(high))
b,a = scipy.signal.iirfilter(order, Wn)
print (a, b) 

wf1 = scipy.signal.lfilter(b, a, wf1)
wf2 = scipy.signal.lfilter(b, a, wf2)
wf3 = scipy.signal.lfilter(b, a, wf3)

numharm = samprate*harmonic_length/harmonic_freq
harm = int(numharm)
han = np.hanning(harm)
for i in range(harm):
	wfm.append(np.sin(2.*np.pi*float(i)*harmonic_freq/samprate))
for j in range(harm):
	wfm[j] = wfm[j]*han[j]

#wf4 = scipy.signal.lfilter(b, a, wf3)
#wfm = scipy.signal.lfilter(b, a, wfm)

# Calculate cross-correlations of master with all three waveforms


cross1 = np.correlate(wf1, wfm,'full')
cross2 = np.correlate(wf2, wfm,'full')
cross3 = np.correlate(wf3, wfm,'full')
#
# Plotting of autocorrelations
#
m0=cross1.max()
t=array('f',[])
for j in range (len(cross1)):
	t.append(float(j)/samprate)
#	cross1[j] = cross1[j]/m0
#	cross2[j] = cross2[j]/m0
#	cross3[j] = cross3[j]/m0

		
if envelope:
	hwf1 = scipy.signal.hilbert(cross1)
	print "Hilbert 1"
	hwf2 = scipy.signal.hilbert(cross2)
	print "Hilbert 2"
	hwf3 = scipy.signal.hilbert(cross3)
	print "Hilbert 3"
#	
	cross1 = np.absolute(hwf1)
	cross2 = np.absolute(hwf2)
	cross3 = np.absolute(hwf3)	

t1 = float(np.argmax(cross1))/float(samprate)
twf=array('f',[])

for j in range (len(wf1)):
	twf.append(float(j)/samprate)	
print "samprate = ", samprate, "length =", len(cross1), "Maximum of CC:", m0


ax1=plt.subplot(711)

plt.plot(t, cross1,t1,1.,'ro')
#ax1.set_ylim([-0.1 ,1.1])
fm = string.Formatter()
label = fm.format(' {0} - Crosscorrelation with {1} Hz {2} cycles',triad, harmonic_freq, harmonic_length)
plt.title(label)

ax1=plt.subplot(712, sharex=ax1)


label = fm.format(' {0} - Crosscorrelation',station2)
t2 = float(np.argmax(cross2))/float(samprate)
#plt.title(label)
plt.plot(t, cross2,t2,1.,'go')
#ax1.set_ylim([-0.1 ,1.1])
fm = string.Formatter()
#plt.title(label)
ax1=plt.subplot(713, sharex=ax1)
t3 = float(np.argmax(cross3))/float(samprate)
plt.plot(t, cross3, t3,1.,'bo')
label = fm.format(' {0} - Crosscorrelation',station3)
#ax1.set_ylim([-0.1 ,1.1])
#plt.title(label)
ax1=plt.subplot(714, sharex=ax1)
w3, = plt.plot(t[0:len(wfm)], wfm)
ax1=plt.subplot(715, sharex=ax1)
w1, = plt.plot(twf, wf1)
ax1=plt.subplot(716, sharex=ax1)
w2, = plt.plot(twf, wf2)
ax1=plt.subplot(717, sharex=ax1)
w3, = plt.plot(twf, wf3)
ylabel = triad+'_'+start_time
figurename = triad+'_'+start_time+'.png'
plt.xlabel(ylabel)
pylab.savefig(figurename)
plt.show()
