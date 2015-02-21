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
import pylab
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

def HydroGroupCCRefine(hag,w1,w2,w3,win_front,win_back,samprate,refine,halfwidth,plot_crossco,plot_detect):
	'''
	Refine the delta-t times obtained from LTA/STA using crosscorrelation 
        and get an updated direction from these.
	'''	
#
# Compute either the crosscorrelation of the signal or the crosscorrelation of the envelopes of the signal in a window set to be win_front seconds before the STA-LTA pick and win_back seconds after.
#
	figc=plt.figure(figsize=(12,6))
#	axc=figc.add_axes([0.1,0.1,0.9,0.9])
	subc1=figc.add_subplot(3,1,1)
	
	label = fm.format('Peak crossco. envelope- {0} detections', len(hag))
	if(refine):
		label = fm.format('Weighted crossco. envelope - {0} detections', len(hag))
	label = triad+' '+label
	plt.title(label)
	subc2=figc.add_subplot(3,1,2, sharex=subc1)
	subc3=figc.add_subplot(3,1,3, sharex=subc1)
	for i in range(len(hag)):		
#
# For each group, get the time window, envelope, and cross-correlation
#
		wf1 = w1[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)]
		wf2 = w2[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)]		
		wf3 = w3[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)]

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
		if(plot_detect):
			figfft = plt.figure(figsize=(10,6))
#			ax = figfft.add_axes([0.1, 0.1, 0.8, 0.8])
			subf1=figfft.add_subplot(1,1,1)						
	                fft1 = abs(scipy.fft(wf1))
			freqs = scipy.fftpack.fftfreq(len(wf1),1./samprate)
			band2050 = 0.
			band5080 = 0.
			nf2050 = 0
			nf5080 = 0
			for f in range(len(freqs)):
				
				if (freqs[f] > 20. and freqs[f] < 50.) :
					nf2050 += 1
					band2050 += fft1[f]
				if (freqs[f] > 50. and freqs[f] < 80.) :
					nf5080 += 1
					band5080 += fft1[f]
			band2050 /= float(nf2050)
			band5080 /= float(nf5080)
			label = fm.format('ratio 5080/2050 is {0}',band5080/band2050)
			plt.title(label)
			subf1.plot(freqs,20.*scipy.log10(fft1),'red')
			

#
# Compute Cross-Correlations
#
		lcor = int(tcor*samprate)
		
		crossco0 = array('f',[])
		crossco1 = array('f',[])
		crossco2 = array('f',[])
		crossco3 = array('f',[])
#		crossco0 = np.correlate(wf1,wf1,'full')
		crossco1 = np.correlate(wf2,wf1,'full')
		#
		#
		#
		print ("Crosscorrelation 1")
		crossco2 = np.correlate(wf3,wf1,'full')
		crossco3 = np.correlate(wf3,wf2,'full')
#
# Hilbert transform
#
		if(envelope):
			hwf1 = scipy.signal.hilbert(crossco1)
			print ("Hilbert 1")
			hwf2 = scipy.signal.hilbert(crossco2)
			print ("Hilbert 2")
			hwf3 = scipy.signal.hilbert(crossco3)
			print ("Hilbert 3")
#
# Envelope
#
			wf1 = np.absolute(hwf1)
			wf2 = np.absolute(hwf2)
			wf3 = np.absolute(hwf3)


		print ("Crosscorrelation 2")
		t0 = float((len(wf1)-1))
		t0 = t0/2.
		t1 = gy.TimeOfMax(wf1)
		t2 = gy.TimeOfMax(wf2)
		t3 = gy.TimeOfMax(wf3)
		dt1 = (t1-t0)/samprate
		dt2 = (t2-t0)/samprate
		dt3 = (t3-t0)/samprate
		ref1=max(wf1)
		ref2=max(wf2)
		ref3=max(wf3)
		time = array('f',[])
		for j in range(len(wf1)):
			time.append(float(j)/samprate - t0/samprate)
		
		
		if (refine):
			maxvalue = wf1.max()
			sumvalue = sum(wf1)/samprate
			ref1=maxvalue/sumvalue
			for j in range (len(wf1)):
				wf1[j] = wf1[j]/sumvalue
			n_it=0
			dtdiff = len(wf1)/samprate
			while (abs(dtdiff) > .5/samprate and n_it < 100):
				t1,std1 = gy.FitLaplace(wf1, samprate, halfwidth, dt1)
				dt1prev=dt1
				dt1 = t1 -dt1 
				dtdiff = dt1-dt1prev
				n_it += 1
			maxvalue = wf2.max()
			sumvalue = sum(wf2)/samprate
			ref2=maxvalue/sumvalue
			for j in range (len(wf2)):
				wf2[j] = wf2[j]/sumvalue	
			n_it=0
			dtdiff = len(wf1)/samprate
			while (abs(dtdiff) > .5/samprate and n_it < 100):
				t2,std2 = gy.FitLaplace(wf2, samprate, halfwidth, dt2)
				dt2prev=dt2
				dt2 = t2 -dt2 
				dtdiff = dt2-dt2prev
				n_it += 1
			maxvalue = wf3.max()
			sumvalue = sum(wf3)/samprate
			ref3=maxvalue/sumvalue
			for j in range (len(wf3)):
				wf3[j] = wf3[j]/sumvalue
			n_it=0
			dtdiff = len(wf1)/samprate
			while (abs(dtdiff) > .5/samprate and n_it < 100):
				t3,std3 = gy.FitLaplace(wf3, samprate, halfwidth, dt3)
				dt3prev=dt3
				dt3 = t3 -dt3 
				dtdiff = dt3-dt3prev
				n_it += 1
		print ("HAG number, HAG:", i, hag[i])
		
		Detnew = hag[i].Det1._replace(dt_cc=hag[i].Det1.time)
		hag[i] = hag[i]._replace(Det1=Detnew)
		Detnew = hag[i].Det2._replace(dt_cc=hag[i].Det2.time+dt1)
		hag[i] = hag[i]._replace(Det2=Detnew)		
		Detnew = hag[i].Det3._replace(dt_cc=hag[i].Det3.time+dt2)
		hag[i] = hag[i]._replace(Det3=Detnew)
		sta1 = hag[i].Det1.sta
		sta2 = hag[i].Det2.sta
		sta3 = hag[i].Det3.sta
		dt = dt1, dt2, dt3, std1, std2, std3
						
		OptDirection, OptVelocity, theta, radii, maxdiff = gy.CompOptDirection(dt,sfile,triad)
		print ("New direction:", OptDirection)
#		print ("Values:", radii)
		hag[i]=hag[i]._replace(direction_cc=OptDirection)
		print ("HAG number, New HAG:", i, hag[i])
		if(plot_crossco):
			
			subc1.plot(time,wf1,'red',dt1,ref1,'bo')
			if (i == 1) :
				label=fm.format('{0}-{1}',sta1,sta2)
				subc1.text(10.,0.7*ref1,label)
				label=fm.format('{0}-{1}',sta1,sta3)
				subc2.text(10.,0.7*ref2,label)
				label=fm.format('{0}-{1}',sta2,sta3)
				subc3.text(10.,0.7*ref3,label)
			subc2.plot(time,wf2,'green',dt2,ref2,'ro')
			subc3.plot(time,wf3,'blue',dt3,ref3,'ro')			

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

snr = config.getfloat('Detect','snr')
t_group = config.getfloat('Detect','time_group')
S_ta = config.getfloat('Detect','STA')
L_ta = config.getfloat('Detect','LTA')
plot_detect = config.getboolean('Detect','plot_detect')
plot_crossco = config.getboolean('Detect','plot_crossco')

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
HydroGroupCCRefine(hag,fltwf1,fltwf2,fltwf3,win_front,win_back,samprate,refine,halfwidth,plot_crossco,plot_detect)

#
#
#

Lshift = int(L_ta*samprate)
Sshift = int(S_ta*samprate)


#
# Plotting of filtered waveforms
#

if(vector_plot == True):
	plt.show()	
	tsta=array('f',[])
	for j in range (len(StaLta1)):
		tsta.append(float(Lshift-Sshift)/samprate+float(j)/samprate)	
	
	ax1=plt.subplot(311)

	pick1 = array('f',[])	
	pick2 = array('f',[])	
	pick3 = array('f',[])

	snr1 = array('f',[])
	nwhale = 0
	for i in range(len(hag)):
		if(hag[i].whale_disc):
			pick1.append(float(Lshift-Sshift)/samprate+hag[i].Det1.time)		
			pick2.append(float(Lshift-Sshift)/samprate+hag[i].Det2.time)		
			pick3.append(float(Lshift-Sshift)/samprate+hag[i].Det3.time)
			snr1.append(snr)
			nwhale += 1
	plt.plot(tsta,StaLta1,'red',pick1,snr1, 'ro')
	plt.plot(tsta,StaLta2,'green',pick2,snr1, 'go')
	plt.plot(tsta,StaLta3,'blue',pick3,snr1, 'bo')

	plt.ylabel("STA/LTA ratio")

	plt.subplot(312,sharex=ax1)
	
	snr2 = array('f',[])
	tflt = array('f',[])
	tpick1 = array('f',[])
	allpick1 = array('f',[])
	for j in range (len(fltwf1)):
		tflt.append(float(j)/samprate)		
	
	for i in range(nwhale):		
		snr2.append(-0.5*max(abs(fltwf1)))

	for i in range(len(det1)):
		allpick1.append(float(Lshift-Sshift)/samprate+det1[i].time)
		tpick1.append(0.)
	for i in range(len(det2)):
		allpick1.append(float(Lshift-Sshift)/samprate+det2[i].time)
		tpick1.append(0.)
	for i in range(len(det3)):
		allpick1.append(float(Lshift-Sshift)/samprate+det3[i].time)
		tpick1.append(0.)
	
	plt.plot(tflt,fltwf1, 'black',allpick1,tpick1, 'ys', pick1, snr2, 'ro')
		
	label = fm.format(' {0} - Filtered',station1)
	plt.ylabel(label)
	ax=plt.subplot(313,sharex=ax1)
        
	X = array('f',[])
	Y = array('f',[])
	u = array('f',[])
	v = array('f',[])
	j=0
	for i in range(len(hag)):
		if(hag[i].whale_disc):
			X.append(pick1[j])
			X.append(pick1[j])
			Y.append(-1.*snr)
			Y.append(snr)
			u.append(math.cos(hag[i].direction))
			v.append(math.sin(hag[i].direction))		
			u.append(math.cos(hag[i].direction_cc))
			v.append(math.sin(hag[i].direction_cc))
			label = fm.format('Angle {0} ',180.*(hag[i].direction-0.5*np.pi)/np.pi)
			ax.text(pick1[j],-0.8*snr,label)
			label = fm.format('Angle {0} ',180.*(hag[i].direction_cc-0.5*np.pi)/np.pi)
			ax.text(pick1[j],1.2*snr,label)
			j += 1

	plt.plot(pick1, snr2, 'ro')
	plt.plot(pick1, snr1, 'bo')
	plt.quiver(X,Y,u,v)
	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)	
	
	plt.ylabel("Incoming direction of signal")
	plt.ylim(-2*snr,2*snr)
	plt.show()
