#!/usr/bin/env Python

# This module reads a configuration file containing 
# the wfdisc tuples, waveform files, site files, and processing parameters.
# The waveform files are written in CSS format four byte integers.
#
# 1 - Reads the waveform and site date
# 2 - Displays the master waveforms 
# 3 - Displays the directional information using either envelope or sta/lta 
#    
# 
# Developped February 2013 - Ronan Le Bras - Gydatos LLC
#

# import all necessary modules

import struct
import math
from   array import *
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from   matplotlib.widgets import Button, Slider
from   collections import namedtuple
import ConfigParser
import string
import time
import calendar
import scipy.signal
import numpy as np
import scipy.io.matlab.byteordercodes
import pylab
import pyaudio
import wave
import glob
import os
import subprocess
import random
import gydatos as gy

Site = namedtuple('Site', ['sta','lat','lon'], verbose=True)

def shiftInPlace(l, n):
	n = n % len(l)
	head = l[:n]
#    l[:n] = []
	l.extend(head)
	shl = l[n:]
	return shl




class Filter:
#
# Callback for the Filter button. 
#
	def filt(self,event):
		w.set_ydata(fltwfm)
		plt.subplot(212,sharex=ax1)
#		Pxx,freqs,bins,im = pylab.specgram(fltwfm,NFFT=512,Fs=samprate,noverlap=500)
		pylab.imshow(filt_im)
		print "In Filter callback!!"
		plt.draw()
	def unfilt(self,event):
		print "In UnFilter callback!!"
		w.set_ydata(wfm)
		plt.subplot(212,sharex=ax1)
#		Pxx,freqs,bins,im = pylab.specgram(wfm,NFFT=512,Fs=samprate,noverlap=500)
		pylab.imshow(im)
		plt.draw()		

# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('Waveforms.cfg')

verbose = config.getint('control_parameters', 'verbose')

wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

start_time = config.get('time_interval','start_time')
end_time = config.get('time_interval','end_time')
time_format = config.get('time_interval','format') 

# station = config.get('time_interval','Station')

triad = config.get('Correlation','Triad')

station1 = triad+'1'
station2 = triad+'2'
station3 = triad+'3'

sta_master   = config.get('Correlation','station_master')
stime_master = config.get('Correlation','start_master')
etime_master = config.get('Correlation','end_master')

order    = config.getint('Filter','order')
low      = config.get('Filter','low')
high     = config.get('Filter','high')
envelope = config.getboolean('Filter','envelope')

snr  = config.getfloat('Detect','snr')
S_ta = config.getfloat('Detect','STA')
L_ta = config.getfloat('Detect','LTA')

sfile    = config.get('Azimuth','sfile')
az_width = config.getfloat('Azimuth','width')
az_type  = config.get('Azimuth','type')
test_case = config.getboolean('Azimuth','test_case')

if verbose > 2:
	print 'Order =', order, 'Low =', low, 'High =', high


# Parse time variables and get start and end epoch time

stime = time.strptime(start_time, time_format)
etime = time.strptime(end_time, time_format)

mstime = time.strptime(stime_master, time_format)
metime = time.strptime(etime_master, time_format)

start_e = calendar.timegm(stime)
end_e   = calendar.timegm(etime)

start_master_e = calendar.timegm(mstime)
end_master_e   = calendar.timegm(metime)


if verbose > 1:
	print "Start time:", stime,'\n'
	print "Epoch start:", start_e
	print "End time:", etime,'\n'
	print "Epoch end:", end_e

samprate, wf1 = gy.ReadWfm(start_e, end_e, station1)
samprate, wf2 = gy.ReadWfm(start_e, end_e, station2)
samprate, wf3 = gy.ReadWfm(start_e, end_e, station3)

#pad=len(wf3)-len(master_wf)
#for i in range (pad):
#	master_wf.append(0.)
# Design digital IIR filter

d1,d2 = gy.ReadSite(sfile, triad)


#
#
#
#
Wn = array('f',[])
Wn.append(float(low))
Wn.append(float(high))
b,a = scipy.signal.iirfilter(order, Wn)
print a, b 

fltwf1 = scipy.signal.lfilter(b, a, wf1)
print "Filtered 1"
fltwf2 = scipy.signal.lfilter(b, a, wf2)
print "Filtered 2"
fltwf3 = scipy.signal.lfilter(b, a, wf3)
print "Filtered 3"
if envelope:
	hwf1 = scipy.signal.hilbert(fltwf1)
	print "Hilbert 1"
	hwf2 = scipy.signal.hilbert(fltwf2)
	print "Hilbert 2"
	hwf3 = scipy.signal.hilbert(fltwf3)
	print "Hilbert 3"
	fltwf1 = np.absolute(hwf1)
	fltwf2 = np.absolute(hwf2)
	fltwf3 = np.absolute(hwf3)
#
# Calculate median level of the envelopes
#
#
sortwf1 = sorted(fltwf1)
medind = int(len(fltwf1)/2)
median1 = sortwf1[medind]
sortwf2 = sorted(fltwf2)
median2 = sortwf2[medind]
sortwf3 = sorted(fltwf3)
median3 = sortwf3[medind]

t1 = [0.,(end_e-start_e)]
med1 = [snr*median1, snr*median1]
med2 = [snr*median2, snr*median2]
med3 = [snr*median3, snr*median3]

StaLta1 = CompStaLta(fltwf1, samprate, S_ta, L_ta)
StaLta2 = CompStaLta(fltwf2, samprate, S_ta, L_ta)
StaLta3 = CompStaLta(fltwf3, samprate, S_ta, L_ta)

if (az_type == 'stalta'):
	t=array('f',[])
	for j in range (len(StaLta1)):
		t.append(float(j)/samprate)		
	print "samprate = ", samprate, "length =", len(fltwf1)
	print "Median1=", median1, "Median2=", median2, "Median3=", median3
	ax1=plt.subplot(311)
	plt.plot(t, StaLta1, 'b')
	fm = string.Formatter()
	label = fm.format('{0} - StaLta',station1)
	plt.title(label)
	plt.subplot(312,sharex=ax1)
	plt.plot(t, StaLta2, 'b')

	label = fm.format('{0} - StaLta',station2)
	plt.title(label)

	plt.subplot(313,sharex=ax1)
	plt.plot(t, StaLta3, 'b')
	plt.ylabel("Correlation coefficient")
	label = fm.format('{0} - StaLta',station3)
	plt.title(label)

if (az_type == 'env'):
	t=array('f',[])
	for j in range (len(fltwf1)):
		t.append(float(j)/samprate)		
	print "samprate = ", samprate, "length =", len(fltwf1)
	print "Median1=", median1, "Median2=", median2, "Median3=", median3
	ax1=plt.subplot(311)
	plt.plot(t, fltwf1, 'b')
	fm = string.Formatter()
	label = fm.format('{0} - StaLta',station1)
	plt.title(label)
	plt.subplot(312,sharex=ax1)
	plt.plot(t, fltwf2, 'b')

	label = fm.format('{0} - StaLta',station2)
	plt.title(label)

	plt.subplot(313,sharex=ax1)
	plt.plot(t, fltwf3, 'b')
	plt.ylabel("Sta Lta")
	label = fm.format('{0} - StaLta',station3)
	plt.title(label)
#
# Computing radial plot values
#

N = 72
theta = np.arange(0.0, 2*np.pi, 2*np.pi/N)
if (az_type == 'env') : 
	print "Computing semblance on envelope\n"
	radii = gy.Semblance(theta, fltwf1,  fltwf2,  fltwf3,  samprate, d1, d2)
elif (az_type == 'stalta'):
	radii = gy.Semblance(theta, StaLta1, StaLta2, StaLta3, samprate, d1, d2)
elif (test_case):
	
	v1 = config.getfloat('Azimuth','theta1')
	dt1,dt2 = CompDt(d1, d2, v1)
	for i in range(len(fltwf1)):
		fltwf1[i] = 0.5*random.random()
		fltwf2[i] = 0.5*random.random()
		fltwf3[i] = 0.5*random.random()
	for i in range(10):
		fltwf1[len(fltwf1)/2-5+i] = 1.
		fltwf2[len(fltwf1)/2-5+i+ int(dt1*samprate)] = 1.
		fltwf3[len(fltwf1)/2-5+i+ int(dt2*samprate)] = 1.
	t=array('f',[])
	for j in range (len(fltwf1)):
		t.append(float(j)/samprate)		
	print "samprate = ", samprate, "length =", len(fltwf1)
	print "Median1=", median1, "Median2=", median2, "Median3=", median3
	ax1=plt.subplot(311)
	plt.plot(t, fltwf1, 'b')
	fm = string.Formatter()
	label = fm.format('{0} - StaLta',station1)
	plt.title(label)
	plt.subplot(312, sharex=ax1)
	plt.plot(t, fltwf2, 'b')

	label = fm.format('{0} - StaLta',station2)
	plt.title(label)

	plt.subplot(313, sharex=ax1)
	plt.plot(t, fltwf3, 'b')
	plt.ylabel("Sta Lta")
	label = fm.format('{0} - StaLta',station3)
	plt.title(label)
#
	radii = gy.Semblance(theta, fltwf1, fltwf2, fltwf3, samprate, d1, d2)

print radii
fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)

width = 2.*np.pi/N
bars = ax.bar(theta, radii, width=width, bottom=0.)
for r,bar in zip(radii, bars):
    bar.set_facecolor( cm.jet(r/10.))
    bar.set_alpha(0.5)

plt.show()

