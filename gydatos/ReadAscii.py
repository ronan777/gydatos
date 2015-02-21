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
import math






class SoundFile:
#
# Play the waveform on audio at 10 times the rate.
# Call back for the button to dump the sound file.
#

	def dumpraw(self,event):
		sound = array('h',[])

		for i in range (len(wfm)):
			sound.append( int (wfm[i] ))
		prefix='raw'+'_'+start_time+'_'+station
		gy.WriteSoundFile(sound,prefix,samprate, 10)

	def dumpfilt(self,event):
		sound = array('h',[])

		for i in range (len(fltwfm)):
			sound.append( int (fltwfm[i] ))
		prefix = 'filt'+'_'+start_time+'_'+station
		gy.WriteSoundFile(sound,prefix,samprate, 10)
		

class Filter:
#
# Callback for the Filter button. 
#
	def filt(self,event):
		w.set_ydata(fltwfm)
		plt.subplot(212,sharex=ax1)

		pylab.imshow(filt_im)
		print ("In Filter callback!!")
		plt.draw()
	def unfilt(self,event):
		print ("In UnFilter callback!!")
		w.set11_ydata(wfm)
		plt.subplot(212,sharex=ax1)

		pylab.imshow(im)
		plt.draw()		

# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('ReadAllFiles_PUBA.cfg')

verbose = config.getint('control_parameters', 'verbose')

wfdisc = config.get('waveforms', 'wfdisc')
file1 = config.get('waveforms', 'ascii_file1')
file2 = config.get('waveforms', 'ascii_file2')

start_time = config.get('time_interval','start_time')
end_time = config.get('time_interval','end_time')
time_format = config.get('time_interval','format') 

station = config.get('time_interval','Station')

order = config.getint('Filter','order')
low = config.get('Filter','low')
high = config.get('Filter','high')

test_case = config.getboolean('Azimuth','test_case')
theta1  = config.getfloat('Azimuth','theta1')
boxcar  = config.getfloat('Azimuth','boxcar')
noise_level = config.getfloat('Azimuth','noise_level')

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


samprate=0.05
t, wfm = gy.Readflat(file1, 31, verbose)
print("File 1. Read samples",len(wfm),"in file", file1)
#
# Master and Main shock
#
wfm1=wfm[0:30000]
t1=t[0:30000]
# CONA
#wfm1=wfm[2*72000+8000:2*72000+10000]
#t1=t[2*72000+8000:2*72000+10000]
tmin=t1[0]
for i in range (len(t1)):
	t1[i]=t1[i]-tmin

t, wf = gy.Readflat(file2, 31, verbose)
print("File 2. Read samples",len(wf),"in file", file2)
wfm2=wf[0:30000]
t2=t[0:30000]
tmin=t2[0]
for i in range (len(t2)):
	t2[i]=t2[i]-tmin
values=file1.split("/")
if(len(values) >1):
	print("Event:",values[3])
	event1 = values[3]
else:
	event1 = values[0]

values=file2.split("/")
if(len(values) >1):
	print("Event:",values[3])
	event2 = values[3]
else:
	event2 = values[0]
dumpfile=open("Waveform_dump",'w')
for i in range (len(wfm1)):
	dumpfile.write("%s\n"% wfm1[i])

dumpfile.close()


#
#
#
if (test_case):

	for i in range (len(wfm2)):
		wfm2[i] = noise_level*np.random.random_sample()
#
# Add master trace in the middle of the trace 
#
	for i in range (len(wfm1)):
		wfm2[len(wfm2)/2-int(len(wfm1)*samprate/2)+i] += wfm1[i]

	for i in range (len(wfm1)):
		wfm2[len(wfm2)/4-int(len(wfm1)*samprate/2)+i] += 0.33*wfm1[i]

	for i in range (len(wfm1)):
		wfm2[3*len(wfm2)/4-int(len(wfm1)*samprate/2)+i] += 0.5*wfm1[i]

# Plotting of waveform and spectrogram

fltwfm1 = scipy.signal.lfilter(b, a, wfm1)
fltwfm2 = scipy.signal.lfilter(b, a, wfm2)

m01=0.
for i in range(len(fltwfm1)):
	m01 += fltwfm1[i]*fltwfm1[i]	
m02=0.
for i in range(len(fltwfm2)):
	m02 += fltwfm2[i]*fltwfm2[i]	

cross1 = np.correlate(fltwfm2, fltwfm1, 'full')
#normcross = gy.equalizeCC(fltwfm2, len(fltwfm1))

#if(len(normcross) != len(cross1)):
#	print("Warning: normcross not the right length", len(cross1), "vs.", len(normcross))

norm=math.sqrt(m01)
t=array('f',[])
for j in range (len(cross1)):
	cross1[j] = cross1[j]/norm
maxcross=cross1.max()

print("Maximum cross-correlation:", maxcross)

ax1 = plt.subplot(211)
#
#
#
w, = plt.plot(t1,wfm1)
fm = string.Formatter()
label = fm.format(' {0} - Raw', event1)
plt.title(label)
ax1=plt.subplot(212, sharex=ax1)
#
#
#
#w, = plt.plot(t1,fltwfm1)
#fm = string.Formatter()
#label = fm.format(' {0} - Filtered - {1}-{2} Nyquist', event1, low, high)
#plt.title(label)
#ax1 = plt.subplot(212, sharex=ax1)
#
#
#
w, = plt.plot(t2,wfm2)
fm = string.Formatter()
label = fm.format(' {0} - Raw', event2)
plt.title(label)
#ax1=plt.subplot(212, sharex=ax1)
#
#
#
#w, = plt.plot(t2,fltwfm2)
#fm = string.Formatter()
#label = fm.format(' {0} - Filtered - {1}-{2} Nyquist', event2, low, high)
#plt.title(label)
#
#
#
#ax1=plt.subplot(413, sharex=ax1)
#tc=array('f',[])
#for i in range(len(cross1)):
#	tc.append(-1.*samprate*len(fltwfm1)+i*samprate)
#w, = plt.plot(tc,cross1)
#fm = string.Formatter()
#label = fm.format('Crosscorrelation {0:.2f}', maxcross)
#plt.title(label)

#ax1=plt.subplot(414, sharex=ax1)
#w, = plt.plot(tc,normcross)
#low_bound=-1.*samprate*len(fltwfm1)
#high_bound=low_bound+samprate*len(cross1)
#ax1.set_xlim([low_bound,high_bound])
#figurename = event1+'_'+event2+'.png'
#pylab.savefig(figurename)
plt.show()
