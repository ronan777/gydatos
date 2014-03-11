#!/usr/bin/env Python

# This module reads a configuration file containing 
# the wfdisc tuples and waveform files in CSS 
# four byte integer s4 big endian format and:
#
# 1 - Displays the content of the waveform 
# 2 - Displays the associated spectrogram on the same x axis as the waveform
# 3 - Dumps a .wav file normalized to an absolute maximum of 30,000 and a ten times
#     accelerated waveform to be able to hear most of the signal.
# 
# Developped February 2013 - copyright(c) Ronan Le Bras - Gydatos LLC
#

# import all necessary modules

import struct
from   array import *
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from   matplotlib.widgets import Button, Slider
import ConfigParser as ConfigParser
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
import pp


def RangeCrossCo(k,w):
	lw = len(w)
	cross = 0.
	for l in range(lw-k):
		cross += w[k+l]*w[l]	
	return cross

def PPAcross(w, length):
	''' 
	Across computes a normalized autocorrelation. The normalization allows to take into account the smaller correlation 
	length at larger lag times
	Input:
		w: input waveform to be correlated with the master waveform mast
		length: length (in number of samples) of the correlation function to be computed.
			The total length will be twice this.
	Output:
		AutoCorrelation function of length 2*length seconds 	
	''' 
	
	across = array('f',[])
	for i in range(length):
		across.append(0.)	

	# Create jobserver
	job_server = pp.Server()

	# Execute the same task with different amount of active workers and measure the time	
		
	start_time = time.time()
	ncpus = job_server.get_ncpus()
	ncpus = ncpus/2
	step = length/ncpus
	
	print ('Starting ', ncpus, ' workers')
	for i in range(step):
		job_server.set_ncpus(ncpus)
		jobs = []
		for k in range(ncpus):

		# Submit a job which will calculate partial autocorrelation 
		# RangeCrossCo - the function
		# (starti, endi) - tuple with arguments for part_sum
		# () - tuple with functions on which function part_sum depends
		# () - tuple with module names which must be imported before part_sum execution

			jobs.append(job_server.submit(RangeCrossCo, (k+i*ncpus, w)))
	   
	        # Retrieve all the results and calculate the resulting array	
		k=0
		for job in jobs:	
			part = job()	
			across[i*ncpus+k] = part
			k=k+1							
				
	        # Print the elapse time
		print ('Time elapsed: ', time.time() - start_time, 's', i, 'of', step)
	
	return across





# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('Waveforms.cfg')

verbose = config.getint('control_parameters', 'verbose')

#wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

start_time = config.get('time_interval','start_time')
end_time = config.get('time_interval','end_time')
time_format = config.get('time_interval','format') 
time_interval = config.getfloat('time_interval','time_interval') 

station = config.get('time_interval','Station')

order = config.getint('Filter','order')
low = config.get('Filter','low')
high = config.get('Filter','high')

S_ta = config.getfloat('Detect','STA')
L_ta = config.getfloat('Detect','LTA')

if verbose > 2:
	print ('Order =', order, 'Low =', low, 'High =', high)

# Design digital IIR filter

Wn = array('f',[])
Wn.append(float(low))
Wn.append(float(high))
b,a = scipy.signal.iirfilter(order, Wn)
print (a, b) 

# Parse time variables and get start and end epoch time

stime = time.strptime(start_time, time_format)
etime = time.strptime(end_time, time_format)

start_epoch = calendar.timegm(stime)
end_epoch   = calendar.timegm(etime)

if verbose > 1:
	print ("Start time:", stime,'\n')
	print ("Epoch start:", start_epoch)
	print ("End time:", etime,'\n')
	print ("Epoch end:", end_epoch)

wfm=array('f',[])

# Get sample rate from wfdisk 

samprate = gy.ReadSrate(wfdisc, station, verbose)

time_interval = max(time_interval, (end_epoch-start_epoch))
for i in range(1+int(time_interval*samprate)):
	wfm.append(0.)
samprate, start_wfm, num_samp = gy.ReadWfm(wfm, wfdisc, start_epoch, end_epoch, station, verbose)

print('Sample rate:',samprate)

t=array('f',[])
for j in range (len(wfm)):
	t.append(float(start_wfm-start_epoch)+float(j))

# Plotting of waveform and cross-correlation

fltwfm = scipy.signal.lfilter(b, a, wfm)
ax1 = plt.subplot(311)


#if(wfm == True):
#	ax1.axis([start_epoch, end_epoch, -1.5*max(wfm),1.5*max(wfm)] ) 
#else:
#	ax1.axis([start_epoch, end_epoch, 0.,1.])

w, = plt.plot(t,fltwfm)
fm = string.Formatter()
label = fm.format(' {0} - Raw',station)
plt.title(label)
plt.subplot(312,sharex=ax1)
StaLta = gy.CompStaLta(L_ta, S_ta, samprate, fltwfm)

print('Length of StaLta, samprate, wfm:', len(StaLta), samprate, len(wfm))

fw, = plt.plot(t[0:len(StaLta)],StaLta[0:len(StaLta)])

label = fm.format('{0} - Filtered',station)
plt.title(label)

plt.subplot(313,sharex=ax1)

#StaLta = gy.CompStaLta(L_ta, S_ta, samprate, fltwfm)

t0=time.time()
print("Time before pylab autocorrelation:",t0)

lags,cross,linecol,b = pylab.acorr(fltwfm, normed=True, maxlags=None, usevlines=True , detrend=mlab.detrend_none)

print("Time elapsed for pylab autocorrelation:",time.time()-t0)

plt.ylabel("Lag in seconds")
#plt.subplot(414,sharex=ax1)
#Pxx,freqs,bins,filt_im = pylab.specgram(fltwfm,NFFT=512,Fs=samprate,noverlap=510)
plt.xlabel(label)
plt.ylabel("Normalized Correlation")

#plt.subplot(414,sharex=ax1)


#mycrossco = PPAcross(fltwfm,len(fltwfm))
#myt = array('f',[])
#for i in range(len(mycrossco)):
#	myt.append(float(i))
#plt.plot(myt, mycrossco)

figurename = station+'_'+start_time+'.png'
pylab.savefig(figurename)
plt.show()
