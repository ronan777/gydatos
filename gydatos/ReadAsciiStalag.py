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
			print(wfm[i])
			sound.append( int (wfm[i] ))
		prefix='raw'+'_'+file1
		gy.WriteSoundFile(sound,prefix,samprate, 5)

	def dumpfilt(self,event):
		sound = array('h',[])

		for i in range (len(fltwfm)):
			sound.append( int (fltwfm[i] ))
		prefix = 'filt'+'_'+file1
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
config.read('Readflat.cfg')

verbose = config.getint('control_parameters', 'verbose')

wfdisc = config.get('waveforms', 'wfdisc')
file1 = config.get('waveforms', 'ascii_file1')

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


samprate=100.
t, wfm = gy.Readflat(file1, 12, verbose)
print("File 1. Read samples",len(wfm),min(wfm),max(wfm))
start_epoch = 0.
start_wfm = 0.
end_epoch   = t[len(t)-1]
print("Limits:",start_wfm, start_epoch, end_epoch)
#
# Master and Main shock
#
#
#
if (test_case):

	for i in range (len(wfm)):
		wfm[i] = noise_level*np.random.random_sample()
#
# Replace traces with boxcar in the middle of the trace 
#
	for i in range(int(boxcar*samprate)):
		wfm[len(wfm)/2-int(boxcar*samprate/2)+i] += 1.

# Plotting of waveform and spectrogram

fltwfm = scipy.signal.lfilter(b, a, wfm)
ax1 = plt.subplot(311)


#if(wfm == True):
#	ax1.axis([start_epoch, end_epoch, -1.5*max(wfm),1.5*max(wfm)] ) 
#else:
#	ax1.axis([start_epoch, end_epoch, 0.,1.])
w, = plt.plot(t,wfm)
fm = string.Formatter()
label = fm.format(' {0} - Raw',station)
plt.title(label)
plt.xlim(float(start_wfm-start_epoch),float(end_epoch-start_epoch))
plt.ylim(min(wfm),max(wfm))
#plt.show()

plt.subplot(312,sharex=ax1)
fw, = plt.plot(t,fltwfm)

label = fm.format('{0} - Filtered',station)
plt.title(label)

plt.subplot(313,sharex=ax1)
Pxx,freqs,bins,im = pylab.specgram(wfm,NFFT=512,Fs=samprate,noverlap=511)
plt.ylabel("Frequency in Hz.")
#plt.subplot(414,sharex=ax1)
#Pxx,freqs,bins,filt_im = pylab.specgram(fltwfm,NFFT=512,Fs=samprate,noverlap=510)
label = fm.format('Seconds after {0}', start_time)
plt.xlabel(label)
plt.ylabel("Frequency in Hz.")
plt.subplots_adjust(bottom=0.1)
callback = Filter()

axfsound  = plt.axes([0.79, 0.025, 0.1, 0.03])
axrsound = plt.axes([0.65, 0.025, 0.1, 0.03])

sound_callback = SoundFile()
bfsound = Button(axfsound, 'Filtered Sound')
bfsound.on_clicked(sound_callback.dumpfilt)
brsound = Button(axrsound, 'Raw Sound')
brsound.on_clicked(sound_callback.dumpraw)

figurename = file1+'.png'
pylab.savefig(figurename)
plt.show()
