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
		w.set_ydata(wfm)
		plt.subplot(212,sharex=ax1)

		pylab.imshow(im)
		plt.draw()		

# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('Readwf.cfg')

verbose = config.getint('control_parameters', 'verbose')

wfdisc = config.get('waveforms', 'wfdisc')

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

wfm=array('f',[])
srate=gy.ReadSrate(wfdisc, station, verbose)
time_interval = end_epoch - start_epoch
for i in range(1+int(time_interval*srate)):
		wfm.append(0.)

if verbose > 1:
	print ("Start time:", stime,'\n')
	print ("Epoch start:", start_epoch)
	print ("End time:", etime,'\n')
	print ("Epoch end:", end_epoch)

samprate, start_wfm, num_samp = gy.ReadWfm(wfm, wfdisc, start_epoch, end_epoch, station, verbose)
print("Read", samprate, start_wfm, num_samp)
dumpfile=open("Waveform_dump",'w')
for i in range (len(wfm)):
	dumpfile.write("%s\n"% wfm[i])

dumpfile.close()
t=array('f',[])
for j in range (len(wfm)):
#	t.append(float(j)/samprate)
	t.append(float(start_wfm)-float(start_epoch)+float(j)/samprate)

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
plt.xlim(float(start_wfm-start_epoch),float(end_epoch-start_epoch))


#if(wfm == True):
#	ax1.axis([start_epoch, end_epoch, -1.5*max(wfm),1.5*max(wfm)] ) 
#else:
#	ax1.axis([start_epoch, end_epoch, 0.,1.])
w, = plt.plot(t,wfm)
fm = string.Formatter()
label = fm.format(' {0} - Raw',station)
plt.title(label)
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

figurename = station+'_'+start_time+'.png'
pylab.savefig(figurename)
plt.show()

