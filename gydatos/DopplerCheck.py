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
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from   matplotlib.widgets import Button, Slider
#import configparser as ConfigParser
import ConfigParser
import string
import time
import calendar
import scipy
import scipy.signal
import scipy.fftpack
import numpy as np
import scipy.io.matlab.byteordercodes
import pylab as pl
#import pyaudio
import wave
import glob
import os
import subprocess
import collections
import gydatos as gy

Detection = collections.namedtuple('Detection',['sta','Det_num','time','dt_cc','duration','freq','whale_disc'], verbose=True)
Hydro_Group = collections.namedtuple('HAG',['triad','HAG_num','whale_disc','Det1','Det2','Det3','direction','direction_cc'], verbose=True)

fm = string.Formatter()

def CheckFreqRatio(signal,samprate,frange,ratio):
	fft1 = abs(scipy.fft(signal))
	freqs = scipy.fftpack.fftfreq(len(signal),1./samprate)
	band1 = 0.
	band2 = 0.
	nf1 = 0
	nf2 = 0
	for f in range(len(freqs)):
				
		if (freqs[f] > frange[0] and freqs[f] < frange[1]) :
			nf1 += 1
			band1 += fft1[f]
		if (freqs[f] > frange[1] and freqs[f] < frange[2]) :
			nf2 += 1
			band2 += fft1[f]
	band1 /= float(nf1)
	band2 /= float(nf2)
	whale = False
	if (band2/band1) > ratio:
		 whale=True
	return whale

def HydroGroupCCRefine(hag,w1,w2,w3,win_front,win_back,samprate,refine,halfwidth,plot_crossco,plot_detect, hanning, plot_waveforms):
	'''
	Refine the delta-t times obtained from LTA/STA using crosscorrelation 
        and get an updated direction from these.
	'''	
#
# Compute either the crosscorrelation of the signal or the crosscorrelation of the envelopes of the signal in a window set to be win_front seconds before the STA-LTA pick and win_back seconds after.
#			
	figfft = plt.figure(figsize=(6,15))
	figamp = plt.figure(figsize=(12,6))
#	ax = figfft.add_axes([0.1, 0.1, 0.8, 0.8])
	subf1=figfft.add_subplot(1,1,1)	
	subf2=figamp.add_subplot(1,1,1)
	const = array('f',[])	
	peak20Hz  = array('f',[])
        amp1_20Hz = array('f',[])
	amp2_20Hz = array('f',[])
	amp3_20Hz = array('f',[])	
	maxamp_20Hz = array('f',[])

	time_20Hz = array('f',[])
	nvalid = 0			
	for i in range(len(hag)):
		peak20Hz.append(0.)
		peak20Hz.append(0.)
		peak20Hz.append(0.)
		
		

	for i in range(len(hag)):		
#
#
#
                Offset = hag[i].Det1.time
	
		wf1 = w1[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)]
		wf2 = w2[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)]		
		wf3 = w3[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)]
#
# Apply a Hanning window
#
		if(hanning):
			han = np.hanning(len(wf1))
			for j in range(len(wf1)):
				wf1[j] = wf1[j]*han[j]
				wf2[j] = wf2[j]*han[j]
				wf3[j] = wf3[j]*han[j]

		frange = 20,50,80
		ratio = 0.2
		whale = False
		whale1 = CheckFreqRatio(wf1, samprate, frange, ratio)
		if (whale1):
			whale = True
		whale2 = CheckFreqRatio(wf2, samprate, frange, ratio)
		if (whale2):
			whale = True
		whale3 = CheckFreqRatio(wf3, samprate, frange, ratio)
		if (whale3):
			whale = True
		print ('Whale discriminant',whale1, whale2, whale3)
		Detnew = hag[i].Det1._replace(whale_disc=whale1)
		hag[i] = hag[i]._replace(Det1=Detnew)
		Detnew = hag[i].Det2._replace(whale_disc=whale2)
		hag[i] = hag[i]._replace(Det2=Detnew)		
		Detnew = hag[i].Det3._replace(whale_disc=whale3)
		hag[i] = hag[i]._replace(Det3=Detnew)
		hag[i] = hag[i]._replace(whale_disc=whale)		
		if(plot_waveforms):			
			figwv = plt.figure(figsize=(8,5))
	#	ax = figfft.add_axes([0.1, 0.1, 0.8, 0.8])
			subf1=figwv.add_subplot(1,1,1)						
			freqs = np.arange(0., len(wf1)/samprate, 1./samprate)	                		
			label = fm.format('HAG number {0}', i)
			plt.title(label)
			print ('Lengths',len(freqs),len(wf1),len(wf2),len(wf3))
			subf1.plot(freqs,wf1,'red',freqs,wf2,'green',freqs,wf3,'blue')
			plt.show()

		if(plot_spectra):	
#			figfft = plt.figure(figsize=(10,6))
#	ax = figfft.add_axes([0.1, 0.1, 0.8, 0.8])
#			subf1=figfft.add_subplot(1,1,1)						
			
			N = int(len(wf1)/2)			
			if(i == 0):
				for c in range(N):
					const.append(0.)

	                fft1 = abs(scipy.fft(wf1))
			freqs = scipy.fftpack.fftfreq(len(wf1),1./samprate)	                
			fft2 = abs(scipy.fft(wf2))	
			fft3 = abs(scipy.fft(wf3))
			label = fm.format('{0} -> {1}', start_time, end_time)	
			constant = 1400.
#
# Check that there is no energy above CONSTANT in interval before and after the expected peak
#	
			nonoise = True
			for c in range(N):
				if(freqs[c] > 15. and freqs[c]  < 18.5):
					if(constant < 200.*scipy.log10(fft1[c])):
						nonoise = False
					if(constant < 200.*scipy.log10(fft2[c])):
						nonoise = False
					if(constant < 200.*scipy.log10(fft3[c])):
						nonoise = False
				if(freqs[c] > 21.5 and freqs[c]  < 25.):
					if(constant < 200.*scipy.log10(fft1[c])):
						nonoise = False
					if(constant < 200.*scipy.log10(fft2[c])):
						nonoise = False
					if(constant < 200.*scipy.log10(fft3[c])):
						nonoise = False
#			if(not nonoise):
#				continue
					
			peak = 0.
			weigh = 0.0000001
			amp1 = 0.
			for c in range(N):
				if(constant >= 200.*scipy.log10(fft1[c])):
					const[c]=200.*scipy.log10(fft1[c])+Offset-constant
					if(freqs[c] > 18. and freqs[c]  < 22.):
						if(amp1 < 20.*scipy.log10(fft1[c])):
							amp1 = 	20.*scipy.log10(fft1[c])
							freq1 = freqs[c]
							
				else:
					const[c]=Offset
					if(freqs[c] > 18. and freqs[c]  < 22.):
						peak += freqs[c]*(200.*scipy.log10(fft1[c])-constant)
						weigh += (200.*scipy.log10(fft1[c])-constant)
						if(amp1 < 20.*scipy.log10(fft1[c])):
							amp1 = 	20.*scipy.log10(fft1[c])
							freq1 = freqs[c]
#			if( weigh == 0.) :
#				continue
			weigh1 = weigh
			peak20Hz[3*i] = freq1
			
			subf1.set_xlim([18. ,22.])
			subf1.set_ylim([0., 14400.])
			subf1.set_title(label)
			print('20Hz peak: r ',peak20Hz[3*i], amp1)
			nvalid += 1
#			subf1.plot(peak20Hz[3*i],Offset+constant-50,'ro')
			subf1.fill_between(freqs[0:N], Offset-constant+200.*scipy.log10(fft1[0:N]), const[0:N], color='red', linewidth=0., edgecolor='white')			
			peak = 0.
			weigh = 0.00001
			amp2 = 0.
			for c in range(N):
				if(constant >= 200.*scipy.log10(fft2[c])):
					const[c]=200.*scipy.log10(fft2[c])+Offset-constant
					if(freqs[c] > 18. and freqs[c]  < 22.):
						if(amp2 < 20.*scipy.log10(fft2[c])):
							amp2 = 	20.*scipy.log10(fft2[c])
							freq2 = freqs[c]
				else:
					const[c]=Offset					
					if(freqs[c] > 18. and freqs[c]  < 22.):
						peak += freqs[c]*(200.*scipy.log10(fft2[c])-constant)
						weigh += (200.*scipy.log10(fft2[c])-constant)
						if(amp2 < 20.*scipy.log10(fft2[c])):
							amp2 = 	20.*scipy.log10(fft2[c])
							freq2 = freqs[c]
#			if( weigh == 0.) :
#				continue			
			peak20Hz[3*i+1] = freq2
			weigh2 = weigh	
				
			print('20Hz peak: g ',peak20Hz[3*i+1], amp2)
			nvalid += 1
#			subf1.plot(peak20Hz[3*i+1],Offset+constant-50.,'go')
			subf1.fill_between(freqs[0:N], Offset-constant+200.*scipy.log10(fft2[0:N]), const[0:N], color='green', linewidth=0., edgecolor='white')			
			peak = 0.
			weigh = 0.00000001
			amp3 = 0.
			for c in range(N):
				if(constant >= 200.*scipy.log10(fft3[c])):
					const[c]=200.*scipy.log10(fft3[c])+Offset-constant
					if(freqs[c] > 18. and freqs[c]  < 22.):
						if(amp3 < 20.*scipy.log10(fft3[c])):
							amp3 = 	20.*scipy.log10(fft3[c])
							freq3 = freqs[c]
				else:
					const[c]=Offset
					if(freqs[c] > 18. and freqs[c]  < 22.):
						peak += freqs[c]*(200.*scipy.log10(fft3[c])-constant)
						weigh += (200.*scipy.log10(fft3[c])-constant)
						if(amp3 < 20.*scipy.log10(fft3[c])):
							amp3 = 	20.*scipy.log10(fft3[c])
							freq3 = freqs[c]
#			if( weigh == 0.) :
#				continue			
			peak20Hz[3*i+2] = freq3
			amp3_20Hz.append(amp3)	
			print('20Hz peak: b',peak20Hz[3*i+2],amp3)	
			time_20Hz.append(Offset)
			amp1_20Hz.append(amp1)
			amp2_20Hz.append(amp2)	
			nvalid += 1	
#			subf1.plot(peak20Hz[3*i+2],Offset+constant-50.,'bo')
			subf1.fill_between(freqs[0:N], Offset-constant+200.*scipy.log10(fft3[0:N]), const[0:N], color='blue', linewidth=0., edgecolor='white')
#			subf1.plot(11.+freqs[0:N],Offset+20.*scipy.log10(fft2[0:N]),'green')
#			subf1.plot(22.+freqs[0:N],Offset+20.*scipy.log10(fft3[0:N]),'blue')
#			plt.show()
		

		if(plot_detect):
			fig = plt.figure(figsize=(6,6))
			ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
			N=len(radii)
			width = 2.*np.pi/N
			bars = ax.bar(theta, radii, width=width, bottom=0.)
			for r,bar in zip(radii, bars):
				bar.set_facecolor( cm.jet(r/10.))
				bar.set_alpha(0.6)
		
			label = fm.format('HAG nb {0} - Refined', i)
			label = triad+' '+label

			plt.title(label)
	plt.grid(True)
#
# Calculate the average peak frequency and plot two lines for the maximum expected Doppler effect.
#
	avg = np.sum(peak20Hz)/float(nvalid)
	print('Average peak: b', avg)
	line1 = array('f',[])
	line1.append(avg-0.2)
	line1.append(avg-0.2)
	line2 = array('f',[])	
	line2.append(avg+0.2)
	line2.append(avg+0.2)
	Y = array('f',[])
	Y.append(0.)
	Y.append(14400.)
	subf1.plot(line1,Y,'k',line2,Y,'k')
	for i in range(len(hag)):	
		maxamp = 0.
		if(amp1_20Hz[i] > maxamp):
			maxamp = amp1_20Hz[i]
			col ='ro'
		if(amp2_20Hz[i] > maxamp):
			maxamp = amp2_20Hz[i]
			col ='go'
		if(amp3_20Hz[i] > maxamp):
			maxamp = amp3_20Hz[i]
			col ='bo'
		subf1.plot(19.,time_20Hz[i],col)
#	subf2.plot(time_20Hz,amp1_20Hz,'ro',time_20Hz,amp2_20Hz,'go',time_20Hz,amp3_20Hz,'bo')				
	plt.show()

def HydroGroupForm(sfile, triad, det1, det2, det3, plot_detect):
	'''
	Form hydroacoustic groups based on detections at the three hydrophones
	'''	
#
# Get maxima of the dts
#
	hag = []
	maxdt1,maxdt2 = gy.CompmaxDt(sfile, triad)
	print ("maxdts",maxdt1,maxdt2)
	hag_num=0

	for i in range(len(det1)):
		for j in range(len(det2)):
			dt1 = det2[j].time-det1[i].time
			if abs(dt1) <= 3.*maxdt1 : 
#
# Look for a third detection in range and declare a group if found
#
				for k in range(len(det3)):
					dt2 = det3[k].time-det1[i].time
					if abs(dt2) <= 3.*maxdt2:
#
# Look for optimum fit direction for these picks
#
						dt3 = det3[k].time-det2[j].time
						std1 = (abs(dt1)+abs(dt1)+abs(dt3))/3.
						std2 = std1 
						std3 = std1
						dt=dt1,dt2, dt3, std1, std2, std3
						
						OptDirection, OptVelocity, theta, radii, maxdiff = gy.CompOptDirection(dt,sfile,triad)
						if(plot_detect):
							fig = plt.figure(figsize=(6,6))
							ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
							N = len(radii)
							width = 2.*np.pi/N
							bars = ax.bar(theta, radii, width=width, bottom=0.)
							for r,bar in zip(radii, bars):
								bar.set_facecolor( cm.jet(r/10.))
								bar.set_alpha(0.5)
	
								
							label = fm.format('HAG nb {0}', hag_num)
							label = triad+' '+label
							plt.title(label)
						whale = False
						newgrp = [triad,hag_num,whale,det1[i],det2[j],det3[k],OptDirection,0.]
						
						print ("Dt vector ", dt)
						print ("New Group", maxdt1, maxdt2,i,j,k,OptDirection,'\n',newgrp,'\n')
						hag.append(Hydro_Group._make(newgrp))
						hag_num+=1
	return hag	
	

def MakeDetect(Det, w, snr, samprate, t_group, sta):
	det = []
	init=0
	det_num=0
	for i in range(1,len(w)):
		t0=float(i)/samprate
		if(w[i] >= snr and w[i] > w[i-1] and init==0):
			init=1
			dw = (snr-w[i-1])/(w[i]-w[i-1])
			if dw != 0. :
				t = t0+1./samprate/dw
			ndet = len(det)
#			print "ndet", ndet, det
			if((ndet==0) or (ndet > 0 and (t-det[ndet-1].time) > t_group)):
				newdet=[sta,det_num,t,0.,0.,0.,False]
				det_num+=1
				print (t0, w[i-1], w[i])
				print (newdet,'\n')
				det.append(Detection._make(newdet))	
		if(w[i] < snr ):
			init=0	
	
	return det		



# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('Waveforms.cfg')

verbose = config.getint('control_parameters', 'verbose')
vector_plot = config.getboolean('control_parameters', 'vector_plot')

#wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

start_time = config.get('time_interval','start_time')
end_time = config.get('time_interval','end_time')
time_format = config.get('time_interval','format') 
time_interval = config.getfloat('time_interval','time_interval') 

# station = config.get('time_interval','Station')

triad = config.get('Correlation','Triad')

station1 = triad+'1'
station2 = triad+'2'
station3 = triad+'3'

correlate = config.getboolean('Correlation','correlation')
sta_master   = config.get('Correlation','station_master')
stime_master = config.get('Correlation','start_master')
etime_master = config.get('Correlation','end_master')
tcor = config.getfloat('Correlation','tcor')
win_front = config.getfloat('Correlation','win_front')
win_back = config.getfloat('Correlation','win_back')

order = config.getint('Filter','order')
low = config.get('Filter','low')
high = config.get('Filter','high')
envelope = config.getboolean('Filter','envelope')
plot_filtered = config.getboolean('Filter','plot_filtered')
hanning = config.getboolean('Filter','hanning')

snr = config.getfloat('Detect','snr')
t_group = config.getfloat('Detect','time_group')
S_ta = config.getfloat('Detect','STA')
L_ta = config.getfloat('Detect','LTA')
plot_detect = config.getboolean('Detect','plot_detect')
plot_crossco = config.getboolean('Detect','plot_crossco')
plot_spectra = config.getboolean('Detect','plot_spectra')
plot_waveforms = config.getboolean('Detect','plot_waveforms')




sfile    = config.get('Azimuth','sfile')
halfwidth = config.getfloat('Azimuth','halfwidth')
az_type  = config.get('Azimuth','type')
test_case = config.getboolean('Azimuth','test_case')
refine = config.getboolean('Azimuth','refine')


if verbose > 2:
	print ('Order =', order, 'Low =', low, 'High =', high)


# Parse time variables and get start and end epoch time

stime = time.strptime(start_time, time_format)
etime = time.strptime(end_time, time_format)

mstime = time.strptime(stime_master, time_format)
metime = time.strptime(etime_master, time_format)

start_e = calendar.timegm(stime)
end_e   = calendar.timegm(etime)

start_master_e = calendar.timegm(mstime)
end_master_e   = calendar.timegm(metime)

end_e = start_e + time_interval

wf1=array('f',[])
wf2=array('f',[])
wf3=array('f',[])
master_wf=array('f',[])

srate=gy.ReadSrate(wfdisc, station1, verbose)
print ("Sampling rate:", srate,'\n')
for i in range(1+int(time_interval*srate)):
		master_wf.append(0.)
samprate, start_master_wf, num_master = gy.ReadWfm(master_wf, wfdisc, start_master_e, end_master_e, sta_master, verbose)

for i in range(1+int(time_interval*samprate)):
		wf1.append(0.)
for i in range(1+int(time_interval*samprate)):
		wf2.append(0.)
for i in range(1+int(time_interval*samprate)):
		wf3.append(0.)		
if verbose > 1:
	print ("Start time:", stime,'\n')
	print ("Epoch start:", start_e)
	print ("End time:", etime,'\n')
	print ("Epoch end:", end_e)

dt1max, dt2max = gy.CompmaxDt(sfile,triad)	

samprate, start_wf1, num_samp1 = gy.ReadWfm(wf1, wfdisc, start_e, end_e, station1, verbose)
samprate, start_wf2, num_samp2 = gy.ReadWfm(wf2, wfdisc, start_e, end_e, station2, verbose)
samprate, start_wf3, num_samp3 = gy.ReadWfm(wf3, wfdisc, start_e, end_e, station3, verbose)
#
#samprate, start_master_wf, num_samp_master = gy.ReadWfm(wfdisc, start_master_e, end_master_e, sta_master, verbose)
#
#
#
Wn = array('f',[])
Wn.append(float(low))
Wn.append(float(high))
b,a = scipy.signal.iirfilter(order, Wn)
print (a, b) 

fltwf1 = scipy.signal.lfilter(b, a, wf1)
print ("Filtered 1")
fltwf2 = scipy.signal.lfilter(b, a, wf2)
print ("Filtered 2")
fltwf3 = scipy.signal.lfilter(b, a, wf3)
print ("Filtered 3")

#
# STA/LTA detector 
#
#
# Detections based on STA/LTA 
#	

StaLta1=gy.CompStaLta(L_ta, S_ta, samprate, fltwf1)
det1=MakeDetect(Detection, StaLta1, snr, samprate, t_group, station1)

StaLta2=gy.CompStaLta(L_ta, S_ta, samprate, fltwf2)
det2=MakeDetect(Detection, StaLta2, snr, samprate, t_group, station2)

StaLta3=gy.CompStaLta(L_ta, S_ta, samprate, fltwf3)
det3=MakeDetect(Detection, StaLta3, snr, samprate, t_group, station3)


hag=HydroGroupForm(sfile, triad, det1, det2, det3, plot_detect)
if(plot_filtered):
	HydroGroupCCRefine(hag,fltwf1,fltwf2,fltwf3,win_front,win_back,samprate,refine,halfwidth,plot_crossco,plot_detect,hanning,plot_waveforms)
else:
	HydroGroupCCRefine(hag,wf1,wf2,wf3,win_front,win_back,samprate,refine,halfwidth,plot_crossco,plot_detect,hanning,plot_waveforms)
#
#
#


