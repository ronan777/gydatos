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
import matplotlib.mpl as mp
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
import random
import collections
import sys

Detection = collections.namedtuple('Detection',
['sta','Det_num','time','dt_cc','duration','freq','whale_disc','pow_15_25','pow_30_50','pow_50_60','spectim'], 
verbose=True)

def CC_snr_atten(snr):
	synmag=array('f',[])
	syndif=array('f',[])
	for i in range(4000):
		synm=i/4000.
		syndif.append(synm)
		rat = math.pow(10.,2.*synm)
		val = 1./math.sqrt(1.+ rat*rat/snr/snr)
		synmag.append(val)	
	return syndif,synmag

def CC_snr_atten_rand(snr):
	synmag=array('f',[])
	syndif=array('f',[])
	for i in range(2000):
		m1 = 4.2*random.random()	
		m2 = 4.2*random.random()
		synm = abs(m1-m2)/abs(m1+m2)
		syndif.append(synm)
#		rat = math.pow(10.,synm)
		a1 = math.pow(10.,m1)
		a2 = math.pow(10.,m2)
		a1a2 = a1*a2
		val =  1./math.sqrt(1.+1./snr/snr/a1/a1)/math.sqrt(1.+1./snr/snr/a2/a2)
#		val = 1./math.sqrt(1.+ rat*rat/snr/snr)
		synmag.append(val)	
	return syndif,synmag


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
config.read('ReadAllFiles.cfg')

verbose = config.getint('control_parameters', 'verbose')
plot_all = config.getboolean('control_parameters','plot_all')
maxcross_thresh = config.getfloat('control_parameters','maxcross_thresh')
mag_thresh = config.getfloat('control_parameters','mag_thresh')

wfdisc = config.get('waveforms', 'wfdisc')
files = config.get('waveforms', 'ascii_files')
magnitudes = config.get('waveforms', 'magnitudes')
pair_plot = config.get('waveforms', 'pair_plot')
pair_plot_corr = config.get('waveforms', 'pair_plot_corr')

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
#
#
#print(files)

file_values=files.split(",")
numfiles=len(file_values)
for i in range(numfiles):
	print(file_values[i])

pair=pair_plot.split(",")
pair_corr=pair_plot_corr.split(",")

mag_values=magnitudes.split(",")
nummag=len(mag_values)
plot_pair=False
if(int(pair[0]) < nummag and int(pair[1]) < nummag):
	plot_pair = True

mags = array('f',[])
mag_diffs = array('f',[])
mag_max = array('f',[])
corr_coefs=array('f',[])
num_clean_mag = 0
for i in range(nummag):
	print(float(mag_values[i]))
	if(mag_values[i] != '-999.'):
		mags.append(float(mag_values[i]))
		num_clean_mag += 1
#
#
#
if(numfiles != num_clean_mag):
	print("ERROR !!! Number of magnitudes does not match number of waveforms",nummag,numfiles,num_clean_mag)
	quit()
#
# Design digital IIR filter
#
Wn = array('f',[])
Wn.append(float(low))
Wn.append(float(high))
b,a = scipy.signal.iirfilter(order, Wn)
print (a, b) 
#
# Parse time variables and get start and end epoch time
#
samprate = 0.01
matrix   = np.zeros((nummag, nummag+2))
mat_dt   = np.zeros((nummag, nummag+2))
colorscale = np.zeros((21,1))
for i in range(20):	
	matrix[i][nummag+1]=i/20.
#
#
filind1=0
filind2=0
c1=array('f',[])
c2=array('f',[])
c3=array('f',[])
c4=array('f',[])
for fil1 in range(nummag):
#
	if(mag_values[fil1] == '-999.'):
		continue
	print('filind1',filind1)
	filind2=0
	for fil2 in range(nummag):
#
#
#	
		print('filind2',filind2)
		if(mag_values[fil2] == '-999.'):
			continue	
				
		file1 = file_values[filind1]
		file2 = file_values[filind2]
		t, wfm = gy.Readflat(file1, 31, verbose)
		print("File 1. Read samples",len(wfm), file1, "first two samples:", wfm[0], wfm[1])
#
#
#wfm1=wfm[3000:13000]
#t1=t[3000:13000]
# CONA
#
#
#wfm1=wfm[50:11900]
#t1=t[50:11900]
# PUBA
		
		wfm1=wfm[0:6200]
		t1=t[0:6200]
		tmin=t1[0]
		for i in range (len(t1)):
			t1[i]=t1[i]-tmin

		t, wf = gy.Readflat(file2, 31, verbose)
		print("File 2. Read samples",len(wf), file2)
		wfm2=wf[0:6200]
		t2=t[0:6200]
		tmin=t2[0]
		for i in range (len(t2)):
			t2[i]=t2[i]-tmin

		values=file1.split("/")
		print("Event:",values[len(values)-1])
		event1 = values[len(values)-1]

		values=file2.split("/")
		print("Event:",values[len(values)-1])
		event2 = values[len(values)-1]

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
#		fltwfm1 = wfm1
#		fltwfm2 = wfm2
		Lta=5.
		Sta=0.5
		snr=4.
		t_group=100.
		station1="XXX"

		StaLta1=gy.CompStaLta(Lta, Sta, 1./samprate, fltwfm1)
		det1=gy.MakeDetect(Detection, StaLta1, snr, 1./samprate, t_group, station1, Lta)

		StaLta2=gy.CompStaLta(Lta, Sta, 1./samprate, fltwfm2)
		det2=gy.MakeDetect(Detection, StaLta2, snr, 1./samprate, t_group, station1, Lta)

		index11=0
		index12=int(25./samprate)
		index21=0	
		index22=int(25./samprate)
		if(len(det1) > 1):
			print("File 1 - NOT A UNIQUE DETECTION !!!") 
			det1[0]._replace(time = 5.)
		if(len(det2) > 1):
			print("File 2 - NOT A UNIQUE DETECTION !!!") 
			det2[0]._replace(time = 5.)

		if(len(det1) > 0 and len(det2) > 0):
			index11 = int((det1[0].time-5.)/samprate)
			index12 = int((det1[0].time+20.)/samprate)
			index21 = int((det2[0].time-5.)/samprate)
			index22 = int((det2[0].time+20.)/samprate)
		m01=0.
		for i in range(len(fltwfm1[index11:index12])):
			m01 += fltwfm1[index11+i]*fltwfm1[index11+i]	
		m02=0.
		for i in range(len(fltwfm2[index21:index22])):
			m02 += fltwfm2[index21+i]*fltwfm2[index21+i]	

		cross1 = np.correlate(fltwfm2[index21:index22], fltwfm1[index11:index12], 'full')
#
#		normcross = gy.equalizeCC(fltwfm2, len(fltwfm1))
#		if(len(normcross) != len(cross1)):
#			print("Warning: normcross not the right length", len(cross1), "vs.", len(normcross))
#
		norm=math.sqrt(m01)*math.sqrt(m02)
		t=array('f',[])
#
#		for j in range (len(cross1)):
#			cross1[j] = cross1[j]/norm/normcross[j]
#
		for j in range (len(cross1)):
			cross1[j] = cross1[j]/norm
		maxcross=cross1.max()
		
		
#		print("Maximum cross-correlation time:", index*samprate)
		corr_coefs.append(maxcross)
		diffmag=abs(mags[filind2]-mags[filind1])
		summag=abs(mags[filind2]+mags[filind1])
		mxmag=max(mags[filind2],mags[filind1])
		if(summag == 0.):
			summag = mxmag
		maxind=np.argmax(cross1)
		dt = -1.*samprate*((len(cross1)-1)/2.)+maxind*samprate
		mag_diffs.append(diffmag/summag)
		mag_max.append(mxmag)
		print("Maximum cross-correlation and index,dtmax", maxcross, dt)
		matrix[fil1][fil2] = maxcross
		tc=array('f',[])
		for i in range(len(cross1)):
			tc.append(-1.*samprate*((len(cross1)-1)/2.)+i*samprate)

		if(abs(dt) < 0.2):
			mat_dt[fil1][fil2] = abs(dt)

		if(plot_pair and int(pair[0]) == fil1 and int(pair[1]) == fil2):
			t3=array('f',[])
			for i in range(len(StaLta1)):
				t3.append(Lta+i*samprate)
			pick1 = array('f',[])	
			pick2 = array('f',[])
			snr1 = array('f',[])
			snr2 = array('f',[])
	
			for i in xrange(len(det1)):
				pick1.append(det1[i].time)
				snr1.append(snr)
			for i in xrange(len(det2)):
				pick2.append(det2[i].time)
				snr2.append(snr)
			ax1 = plt.subplot(411)
			plt.plot(t1,fltwfm1, pick1,snr1,'ro')
			plt.subplot(412, sharex=ax1)
			plt.plot(t2,fltwfm2,pick2,snr2,'ro')
			plt.subplot(413, sharex=ax1)
			delta = int(Lta/samprate)
			plt.plot(t3,StaLta1, pick1,snr1,'ro')
			plt.subplot(414, sharex=ax1)
			plt.plot(t3,StaLta2, pick2,snr2,'ro')
			fm = string.Formatter()
			label = fm.format('Aftershocks {0} (M={1:.1f}) and {2} (M={3:.1f}) - CC={4:.2f}', fil1, mags[filind1],fil2, mags[filind2],maxcross)
			plt.title(label)
			plt.show()

		if(int(pair_corr[0]) == fil1 and int(pair_corr[1]) == fil2):
			ax1 = plt.subplot(111)
			plt.plot(tc,cross1)
			fm = string.Formatter()
			label = fm.format('Aftershocks {0} (M={1:.1f}) and {2} (M={3:.1f}) - CC={4:.2f}', fil1, mags[filind1],fil2, mags[filind2],maxcross)
			plt.title(label)
			plt.show()
	#
	#
	#
	#
	#

		if(plot_all and maxcross > maxcross_thresh and mags[fil1] > mag_thresh and mags[fil2] > mag_thresh):

	#
	#
			ax1=plt.subplot(111)
			
			w, = plt.plot(tc,cross1)
			fm = string.Formatter()
			label = fm.format('Crosscorrelation {0:.2f}', maxcross)
			plt.title(label)
			plt.ylim([-1.0 ,1.1])
			plt.xlim([-1. ,1.])

		filind2 += 1
	filind1 += 1
plt.show()
plt.figure()
#print('pair[0]=',pair[0])
#for i in range(nummag):
#	c1.append(matrix[int(pair[0])][i])
#	c2.append(matrix[int(pair[1])][i])
#	c3.append(float(i)/float(nummag))
#	c4.append(float(i)/float(nummag))
#print('Correlation coefficienst for two largest earthquakes:',c1,c2)				
#plt.plot(c1,c2,'ro',c3,c4,'bo')
#plt.show()

cmap_gr = mp.colors.ListedColormap([[0.,0.,0.],[0.2,0.2,0.2],[0.2,0.2,0.3],[.2,0.3,.3],[0.3,0.3,0.3],[.3,0.3,0.4],[1.,1.,.0],[1.,.0,.0]])
#matrix.sort()
print(matrix)
plt.figure()
#imgplot = plt.imshow(matrix, cmap=cmap_gr) 
imgplot = plt.imshow(matrix, cmap='hot') 
imgplot.set_yticklabels=(['0.','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'])
imgplot.set_interpolation('nearest')
#imgplot.set_cmap('hot')
plt.show()
#sclplot= plt.imshow(colorscale)
#imgplot.set_interpolation('nearest')
#imgplot.set_cmap('hot')
#plt.colorbar()
#plt.show()
print(mat_dt)
imgplot = plt.imshow(mat_dt, cmap=cmap_gr)
imgplot.set_interpolation('nearest')
imgplot.set_cmap('hot')
plt.show()
#syndif5,synmag5=CC_snr_atten_rand(5.)
#syndif3,synmag3=CC_snr_atten_rand(3.)
#syndif1,synmag1=CC_snr_atten_rand(1.)
syndif5,synmag5=CC_snr_atten_rand(.5)
#syndif25,synmag25=CC_snr_atten_rand(.25)

plt.plot(mag_diffs,corr_coefs,'ro',syndif5,synmag5,'b.')
plt.xlim([-0.2, 1.2])
plt.ylim([0. ,1.1])
plt.xlabel("Normalized Magnitude Differences")
plt.ylabel("Zero-lag Cross-Correlation")
plt.show()
plt.plot(mag_max,corr_coefs,'bo')
plt.xlim([-0.2, 4.5])
plt.ylim([0. ,1.1])
plt.xlabel("Maximum Magnitude")
plt.ylabel("Zero-lag Cross-Correlation")
plt.show()


