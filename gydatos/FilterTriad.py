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
from   matplotlib.widgets import Button, Slider
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
import gydatos as gy


def CorrMaster(w, mast, length):

	lmast = len(mast)
	lw = len(w)
	crossco = array('f',[])
	
# Compute the energy in master waveform

	en_master = 0.
	for i in range (lmast):
		en_master += mast[i]*mast[i]	
#		
# Compute the energy in waveform 
#
	en_master = math.sqrt(en_master)
	en_w = 0.
	for j in range (len(w)):
		en_w += w[j]*w[j]
	en_w = math.sqrt(en_w)	
#		
# Compute the energy in waveform 
#	
	cross = 0.
	for k in range (length):
		if (k%10 == False): 
			print k
		for l in range((len(w)-k)):
			cross += w[k+l]*mast[l]
		cross = cross/en_w/en_master
		crossco.append(cross)
	array.reverse(crossco)
	for k in range (length):
		if (k%10 == False):
			print k
		for l in range(1,(len(w)-k)):
			cross += mast[k+l]*w[l]
		cross = cross/en_w/en_master
		crossco.append(cross)
	array.reverse(crossco)
	print crossco
	return crossco

def TimeOfMax(w):
	maximum=0.
	maxindex=0
	for i in range(len(w)):
		if(abs(w[i]) > maximum):
			maximum=abs(w[i])
			maxindex=i
	return maxindex



class SoundFile:
#
# Play the waveform on audio at 10 times the rate.
# Call back for the button to dump the sound file.
#

	def dumpraw(self,event):
		sound = array('h',[])

		for i in range (len(wfm)):
			sound.append( int (wfm[i] ))
		gy.WriteSoundFile(sound,'raw')

	def dumpfilt(self,event):
		sound = array('h',[])

		for i in range (len(fltwfm)):
			sound.append( int (fltwfm[i] ))
		gy.WriteSoundFile(sound,'filt')
		

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
tcor = config.getfloat('Correlation','tcor')

order = config.getint('Filter','order')
low = config.get('Filter','low')
high = config.get('Filter','high')
envelope = config.getboolean('Filter','envelope')

snr = config.getfloat('Detect','snr')
S_ta = config.getfloat('Detect','STA')
L_ta = config.getfloat('Detect','LTA')

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
samprate, master_wf = gy.ReadWfm(start_master_e, end_master_e, sta_master)
#pad=len(wf3)-len(master_wf)
#for i in range (pad):
#	master_wf.append(0.)
# Design digital IIR filter

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

tmed = [0.,(end_e-start_e)]
med1 = [snr*median1, snr*median1]
med2 = [snr*median2, snr*median2]
med3 = [snr*median3, snr*median3]
#
#
# Find peak correlation of envelopes
#
#
lcor=int(tcor*samprate)
crossco0 = array('f',[])
crossco1 = array('f',[])
crossco2 = array('f',[])
crossco3 = array('f',[])
crossco0 = gy.CorrMaster(fltwf1,fltwf1,lcor)
crossco1 = gy.CorrMaster(fltwf1,fltwf2,lcor)
#
#
#
print "Crosscorelation 1"
crossco2 = gy.CorrMaster(fltwf1,fltwf3,lcor)
print "Crosscorelation 2"
t0 = gy.TimeOfMax(crossco0)
t1 = gy.TimeOfMax(crossco1)
t2 = gy.TimeOfMax(crossco2)
dt1 = (t1-t0)/samprate
dt2 = (t2-t0)/samprate
crossco3 = gy.CorrMaster(fltwf2,fltwf3,lcor)
print "Crosscorelation 3"
t3 = gy.TimeOfMax(crossco3)
dt3= (t3-t0)/samprate
diffdt = dt3-(dt2+dt1)
print "t0, t1, t2, t3", t0, t1, t2, t3
print "dt1, dt2, dt3, diffdt", dt1, dt2, dt3, diffdt

#
# STA/LTA detector 
#
# Initiate
#

#
# Plotting of filtered waveforms
#
t=array('f',[])
for j in range (len(fltwf1)):
	t.append(float(j)/samprate)		
print "samprate = ", samprate, "length =", len(fltwf1)
print "Median1=", median1, "Median2=", median2, "Median3=", median3
ax1=plt.subplot(331)
plt.plot(t, fltwf1, 'b', tmed, med1,'r')
fm = string.Formatter()
label = fm.format(' {0} - Filtered',station1)
plt.title(label)
plt.subplot(334,sharex=ax1)
plt.plot(t, fltwf2, 'b', tmed, med2,'r')

label = fm.format('{0} - Filtered',station2)
plt.title(label)

plt.subplot(337,sharex=ax1)
plt.plot(t, fltwf3, 'b',tmed, med3,'r')
plt.ylabel("Correlation coefficient")
label = fm.format('{0} - Filtered',station3)
plt.title(label)
plt.subplot(332,sharex=ax1)
StaLta1=gy.CompStaLta(L_ta, S_ta, samprate, fltwf1)

Lshift = int(L_ta*samprate)
Sshift = int(S_ta*samprate)

tsta=array('f',[])
for j in range (len(StaLta1)):
	tsta.append(float(Lshift-Sshift)/samprate+float(j)/samprate)	
plt.plot(tsta,StaLta1)
label = fm.format('Seconds after {0}',start_time)
plt.xlabel(label)
plt.ylabel("StaLta ratio")

plt.subplot(335,sharex=ax1)
StaLta2=gy.CompStaLta(L_ta, S_ta, samprate, fltwf2)
	
plt.plot(tsta,StaLta2)
label = fm.format('Seconds after {0}',start_time)
plt.xlabel(label)
plt.ylabel("StaLta ratio")

plt.subplot(338,sharex=ax1)
StaLta3=gy.CompStaLta(L_ta, S_ta, samprate, fltwf3)
	
plt.plot(tsta,StaLta3)
label = fm.format('Seconds after {0}',start_time)
plt.xlabel(label)
plt.ylabel("StaLta ratio")
plt.subplot(333,sharex=ax1)
tcor=array('f',[])
for j in range (len(crossco1)):
	tcor.append(-0.5*float(len(crossco1)-1)/samprate+float(j)/samprate)	
	
plt.plot(tcor,crossco1)
label = fm.format('Seconds after {0}',start_time)
plt.xlabel(label)
plt.ylabel("StaLta ratio")
plt.subplot(336,sharex=ax1)
	
plt.plot(tcor,crossco2)
label = fm.format('Seconds after {0}',start_time)
plt.xlabel(label)
plt.ylabel("StaLta ratio")
plt.subplot(339,sharex=ax1)
	
plt.plot(tcor,crossco3)
label = fm.format('Seconds after {0}',start_time)
plt.xlabel(label)
plt.ylabel("StaLta ratio")

#
#plt.subplots_adjust(bottom=0.2)
#
figurename = triad+'_'+start_time+'.png'
pylab.savefig(figurename)
plt.show()
