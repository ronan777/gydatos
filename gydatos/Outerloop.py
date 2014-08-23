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
import matplotlib.mpl as mp
import matplotlib.colorbar as col
import matplotlib.cm as cm
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
import pywt as wlt

Detection = collections.namedtuple('Detection',
['sta','Det_num','time','dt_cc','duration','freq','snr','whale_disc','reclev','pow_15_25','pow_30_50','pow_50_60','spectim'], 
verbose=True)
Hydro_Group = collections.namedtuple('HAG',['triad','HAG_num','whale_disc','Det1','Det2','Det3','direction','absfit','direction_cc','relfit','direction_stk','stkfit'], 
verbose=True)

fm = string.Formatter()


NAN=999.

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
						
						OptDirection,OptVelocity,theta,radii,maxdiff = gy.CompOptDirection(dt,sfile,triad)
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
						newgrp = [triad,hag_num,whale,det1[i],det2[j],det3[k],OptDirection,maxdiff,NAN,0.,NAN,0.]
						
						print ("Dt vector ", dt)
						print ("New Group", maxdt1, maxdt2,i,j,k,OptDirection,'\n',newgrp,'\n')
						hag.append(Hydro_Group._make(newgrp))
						hag_num+=1
	return hag	
	

#def MakeDetect(Det, w, snr, samprate, t_group, sta, shift):

#	det = []
#	spectim=[]
#	init=0
#	det_num=0
#	for i in range(1,len(w)):
#		t0=float(i)/samprate
#		if(w[i] >= snr and w[i] > w[i-1] and init==0):
#			init=1
#			dw = (snr-w[i-1])/(w[i]-w[i-1])
#			if dw != 0. :
#				t = t0+shift+dw/samprate
#			ndet = len(det)
#			print "ndet", ndet, det
#			if((ndet==0) or (ndet > 0 and (t-det[ndet-1].time) > t_group)):
#				newdet=[sta,det_num,t,0.,0.,0.,False,0.,0.,0.,spectim]
#				det_num+=1
#				print (t0, w[i-1], w[i])
#				print (newdet,'\n')
#				det.append(Detection._make(newdet))	
#		if(w[i] < snr ):
#			init=0	
#	
#	return det		



# Read configuration parameters

level1_psi=False
level2_psi=False
level3_psi=False
level4_psi=False
level5_psi=False
level6_psi=False
level7_psi=False
level1_phi=False
level2_phi=False
level3_phi=False
level4_phi=False
level5_phi=False
level6_phi=False
level7_phi=False

config = ConfigParser.RawConfigParser()
config.read('Waveforms.cfg')

verbose = config.getint('control_parameters', 'verbose')
vector_plot = config.getboolean('control_parameters', 'vector_plot')

#wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

start_time  = config.get('time_interval','start_time')
end_time    = config.get('time_interval','end_time')
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
whole_series_spectrum = config.getboolean('Filter','whole_series_spectrum')
band_pass = config.getboolean('Filter','band_pass')
wavelet_transform = config.getboolean('Filter','wavelet_transform')
level_psi = config.getboolean('Filter','level_psi')
level_phi = config.getboolean('Filter','level_phi')

wlt_type = config.get('Filter','wavelet_type')
wlt_side = config.get('Filter','wavelet_side')
add_levels = config.get('Filter','wavelet_add')
start_level = config.getint('Filter','start_level')

snr = config.getfloat('Detect','snr')
t_group = config.getfloat('Detect','time_group')
S_ta = config.getfloat('Detect','STA')
L_ta = config.getfloat('Detect','LTA')
plot_detect = config.getboolean('Detect','plot_detect')
plot_waveform = config.getboolean('Detect','plot_waveform')
plot_crossco = config.getboolean('Detect','plot_crossco')
stack = config.getboolean('Detect','stack')
plot_stack = config.getboolean('Detect','plot_stack')
plot_spect = config.getboolean('Detect','plot_spect')
plot_hist = config.getboolean('Detect','plot_hist')
plot_spectim = config.getboolean('Detect','plot_spectim')
plot_power = config.getboolean('Detect','plot_power')
plot_ratios = config.getboolean('Detect','plot_ratios')


sfile    = config.get('Azimuth','sfile')
halfwidth = config.getfloat('Azimuth','halfwidth')
az_type  = config.get('Azimuth','type')
test_case = config.getboolean('Azimuth','test_case')
theta1  = config.getfloat('Azimuth','theta1')
boxcar  = config.getfloat('Azimuth','boxcar')
noise_level = config.getfloat('Azimuth','noise_level')

refine_crossco = config.getboolean('Azimuth','refine_crossco')
percent = config.getfloat('Azimuth','percent_top')
minfit = config.getfloat('Azimuth','thresh_fit')


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



if verbose > 1:
	print ("Start time:", stime,'\n')
	print ("Epoch start:", start_e)
	print ("End time:", etime,'\n')
	print ("Epoch end:", end_e)

master_wf = array('f',[])
wf1=array('f',[])
wf2=array('f',[])
wf3=array('f',[])
add1 = array('f',[])
add2 = array('f',[])
add3 = array('f',[])
dt1max, dt2max = gy.CompmaxDt(sfile,triad)

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
#
# Loop over interval defined by variable time_interval
#
#

num_intervals = int(math.floor((end_e-start_e)/time_interval))
	
for time_int in xrange(num_intervals):

	label = fm.format('{0}', time_int)	



	end_e = start_e + time_interval
	lab_time=time.gmtime(start_e)
	pr_time = time.strftime(time_format, lab_time)
	print(pr_time)
	figurename = triad+'_'+pr_time+'.png'
	print(figurename)

	samprate, start_wf1, num_samp1 = gy.ReadWfm(wf1, wfdisc, start_e, end_e, station1, verbose)
	samprate, start_wf2, num_samp2 = gy.ReadWfm(wf2, wfdisc, start_e, end_e, station2, verbose)
	samprate, start_wf3, num_samp3 = gy.ReadWfm(wf3, wfdisc, start_e, end_e, station3, verbose)

	
	#
	#
	#
	start_e += time_interval
	if (test_case):
		d1,d2 = gy.ReadSite(sfile, triad, 0)
		theta1 = (theta1+90.)*np.pi/180.
		dt1,dt2 = gy.CompDt(d1,d2,theta1,1.48)
		for i in xrange (len(wf1)):
			wf1[i] = noise_level*np.random.random_sample()
			wf2[i] = noise_level*np.random.random_sample()
			wf3[i] = noise_level*np.random.random_sample()
	#
	# Replace traces with boxcar in the middle of the trace 
	#

		for i in xrange(int(boxcar*samprate)):
			wf1[len(wf1)/2-int(boxcar*samprate/2)+i] += 1.
			wf2[len(wf2)/2-int(boxcar*samprate/2)+i+int(dt1*samprate)] += 1.
			wf3[len(wf3)/2-int(boxcar*samprate/2)+i+int(dt2*samprate)] += 1.

	if(whole_series_spectrum):
		figfft = plt.figure(figsize=(10,6))
		subf1 = figfft.add_subplot(1,1,1)	
	        fft1 = abs(scipy.fft(wf1))/samprate
		freqs = scipy.fftpack.fftfreq(len(wf1),1./samprate)		
		label = 'Spectrum '+pr_time
		print(label)
		plt.title(label)
		ind0=len(wf1)/2
		print(ind0,freqs[0])
		subf1.plot(freqs[0:len(wf1)/2],20.*scipy.log10(fft1[0:len(wf1)/2]),'red')
		plt.show()
	if(band_pass):
		Wn = array('f',[])
		Wn.append(float(low))
		Wn.append(float(high))
		b,a = scipy.signal.iirfilter(order, Wn)
		print (a, b) 
	#	print("Expected sample length:", int(1+(end_e-start_e)*samprate))
		print("Read", num_samp1, num_samp2, num_samp3, "samples")
		fltwf1 = scipy.signal.lfilter(b, a, wf1)
		print ("Filtered 1, length:", len(wf1))
		fltwf2 = scipy.signal.lfilter(b, a, wf2)
		print ("Filtered 2, length:", len(wf2))
		fltwf3 = scipy.signal.lfilter(b, a, wf3)
		print ("Filtered 3, length:", len(wf3))
		n1 = len(wf1)
	elif(wavelet_transform):
		w = wlt.Wavelet(wlt_type)
		N=1
		n=1
		while (len(wf1) > N):
			N = N*2
			if(len(wf1) < N):
				for i in range(N-len(wf1)):
					wf1.append(0.)
					wf2.append(0.)
					wf3.append(0.)
			len1 = wlt.swt_max_level(N)
			n1 = len(wf1)
			len3 = wlt.swt_max_level(n1)
			len2 = wlt.dwt_max_level(N, w.dec_len)
			print(N, 'swt:',len1,'dwt:',len2,len3,n1)
			n += 1		
		for i in range(N):
			add1.append(0.)
			add2.append(0.)
			add3.append(0.)

		[(cA7, cD7),(cA6, cD6),(cA5, cD5),(cA4, cD4),(cA3, cD3),(cA2, cD2), (cA1, cD1)] = wlt.swt(wf1, wlt_type, 7, start_level=start_level)
		if(len(add_levels) > 0):
#			figadd = plt.figure(figsize=(10,6))
#			ax1=plt.subplot(511)
#			plt.title(label)
			for level in range(len(add_levels)):
				print("Level",add_levels[level])
				if(int(add_levels[level]) == 0):
					print("Adding cD7")		
					for i in range(N):
						add1[i] += cD7[i]
				if(int(add_levels[level]) == 1):
					print("Adding cD6")		
					for i in range(N):
						add1[i] += cD6[i]
				if(int(add_levels[level]) == 2):
					print("Adding cD5")		
					for i in range(N):
						add1[i] += cD5[i]
				if(int(add_levels[level]) == 3):
					print("Adding cD4")		
					for i in range(N):
						add1[i] += cD4[i]
				if(int(add_levels[level]) == 4):
					print("Adding cD3")		
					for i in range(N):
						add1[i] += cD3[i]
				if(int(add_levels[level]) == 5):
					print("Adding cD2")		
					for i in range(N):
						add1[i] += cD2[i]
				if(int(add_levels[level]) == 6):
					print("Adding cD1")		
					for i in range(N):
						add1[i] += cD1[i]
#			plt.plot(tt,add1)
			fltwf1 = add1
#			plt.show()		
		
		[(cA7, cD7),(cA6, cD6),(cA5, cD5),(cA4, cD4),(cA3, cD3),(cA2, cD2), (cA1, cD1)] = wlt.swt(wf2, wlt_type, 7, start_level=start_level)
		if(len(add_levels) > 0):
#			figadd = plt.figure(figsize=(10,6))
#			ax1=plt.subplot(412, sharex=ax1)
#			plt.title(label)
			for level in range(len(add_levels)):
				print("Level",add_levels[level])
				if(int(add_levels[level]) == 0):
					print("Adding cD7")		
					for i in range(N):
						add2[i] += cD7[i]
				if(int(add_levels[level]) == 1):
					print("Adding cD6")		
					for i in range(N):
						add2[i] += cD6[i]
				if(int(add_levels[level]) == 2):
					print("Adding cD5")		
					for i in range(N):
						add2[i] += cD5[i]
				if(int(add_levels[level]) == 3):
					print("Adding cD4")		
					for i in range(N):
						add2[i] += cD4[i]
				if(int(add_levels[level]) == 4):
					print("Adding cD3")		
					for i in range(N):
						add2[i] += cD3[i]
				if(int(add_levels[level]) == 5):
					print("Adding cD2")		
					for i in range(N):
						add2[i] += cD2[i]
				if(int(add_levels[level]) == 6):
					print("Adding cD1")		
					for i in range(N):
						add2[i] += cD1[i]
#			plt.plot(tt,add2)
			fltwf2 = add2
		[(cA7, cD7),(cA6, cD6),(cA5, cD5),(cA4, cD4),(cA3, cD3),(cA2, cD2), (cA1, cD1)] = wlt.swt(wf3, wlt_type, 7, start_level=start_level)
		if(len(add_levels) > 0):
#		figadd = plt.figure(figsize=(10,6))
#		ax1=plt.subplot(413, sharex=ax1)
#		plt.title(label)
			for level in range(len(add_levels)):
				print("Level",add_levels[level])
				if(int(add_levels[level]) == 0):
					print("Adding cD7")		
					for i in range(N):
						add3[i] += cD7[i]
				if(int(add_levels[level]) == 1):
					print("Adding cD6")		
					for i in range(N):
						add3[i] += cD6[i]
				if(int(add_levels[level]) == 2):
					print("Adding cD5")		
					for i in range(N):
						add3[i] += cD5[i]
				if(int(add_levels[level]) == 3):
					print("Adding cD4")		
					for i in range(N):
						add3[i] += cD4[i]
				if(int(add_levels[level]) == 4):
					print("Adding cD3")		
					for i in range(N):
						add3[i] += cD3[i]
				if(int(add_levels[level]) == 5):
					print("Adding cD2")		
					for i in range(N):
						add3[i] += cD2[i]
				if(int(add_levels[level]) == 6):
					print("Adding cD1")		
					for i in range(N):
						add3[i] += cD1[i]
#				plt.plot(tt,add3)
				fltwf3 = add3
#
# STA/LTA detector 
#
#
# Detections based on STA/LTA 
#	

	StaLta1 = gy.CompStaLta(L_ta, S_ta, samprate, fltwf1[0:n1])
	det1 = gy.MakeDetect(Detection, StaLta1, fltwf1, snr, samprate, t_group, station1, L_ta)

	StaLta2 = gy.CompStaLta(L_ta, S_ta, samprate, fltwf2[0:n1])
	det2 = gy.MakeDetect(Detection, StaLta2, fltwf2, snr, samprate, t_group, station2, L_ta)

	StaLta3 = gy.CompStaLta(L_ta, S_ta, samprate, fltwf3[0:n1])
	det3 = gy.MakeDetect(Detection, StaLta3, fltwf3, snr, samprate, t_group, station3, L_ta)
	
	hag = HydroGroupForm(sfile, triad, det1, det2, det3, plot_detect)

#	gy.ReceivedLevel(fltwf1, det1, samprate, t_group)
#	gy.ReceivedLevel(fltwf2, det2, samprate, t_group)
#	gy.ReceivedLevel(fltwf3, det3, samprate, t_group)

	gy.WhaleDiscriminant(hag,fltwf1,fltwf2,fltwf3,sfile,win_front,win_back,samprate)

	if(refine_crossco):
		gy.HydroGroupCCRefine(hag, fltwf1, fltwf2, fltwf3, sfile, triad, win_front, win_back,
			samprate, tcor, envelope, refine_crossco, halfwidth,
			plot_waveform, plot_crossco, plot_detect, plot_spect)
	if(stack):
		gy.HydroGroupSTRefine(hag, fltwf1, fltwf2, fltwf3, sfile, triad, win_front, win_back,
			samprate, float(low)*samprate*.5, float(high)*samprate*.5,
			plot_stack, plot_spectim, percent, test_case)

#
#
#

	Lshift = int(L_ta*samprate)
	Sshift = int(S_ta*samprate)

#
# Plotting of filtered and processed waveforms
#

	tsta=array('f',[])
	tsta2=array('f',[])
	tsta3=array('f',[])
	for j in xrange (len(StaLta1)):
		tsta.append(float(Lshift)/samprate+float(j)/samprate)	
	
	f = plt.figure(figsize=(11,8))
	ax1=plt.subplot(512)

	pick1 = array('f',[])	
	pick2 = array('f',[])	
	pick3 = array('f',[])

	snr1 = array('f',[])
	nwhale = 0
	for i in xrange(len(hag)):
#		if(hag[i].whale_disc):
			pick1.append(hag[i].Det1.time)		
			pick2.append(hag[i].Det2.time)		
			pick3.append(hag[i].Det3.time)
			snr1.append(snr)
			nwhale += 1
	plt.plot(tsta,StaLta1,'red',  pick1,snr1, 'ro')

#	print(len(tsta),len(StaLta1),len(StaLta2),len(StaLta3))	

	plt.plot(tsta,StaLta2,'green',pick2,snr1, 'go')
	plt.plot(tsta,StaLta3,'blue' ,pick3,snr1, 'bo')

	plt.xlabel("STA/LTA ratio")

	plt.subplot(511,sharex=ax1)
	
	snr2 = array('f',[])
	snr3 = array('f',[])
	snr1 = array('f',[])
	tflt = array('f',[])
	tpick1 = array('f',[])
	allpick1 = array('f',[])
	for j in xrange (len(fltwf1)):
		tflt.append(float(j)/samprate)		
	hmax = 0.5*max(fltwf1)
	
	for i in xrange(len(pick1)):
		snr1.append(-1.*hmax)
	for i in xrange(len(pick1)):
		snr2.append(0.)
	for i in xrange(len(pick1)):
		snr3.append(hmax)
		
	
	plt.plot(tflt, fltwf1, 'black', pick1, snr1, 'ro', pick2, snr2, 'ys', pick3, snr3, 'ys')
		
	label = fm.format(' {0} - Filtered',station1)
	plt.xlabel(label)
	ax=plt.subplot(513,sharex=ax1)
        
	X = array('f',[])
	Y = array('f',[])
	u = array('f',[])
	v = array('f',[])
	C = array('c',[])
	j = 0

	for i in xrange(len(hag)):
#		if(hag[i].whale_disc and (hag[i].absfit > minfit or hag[i].relfit > minfit)):
		if(hag[i].absfit > minfit or hag[i].relfit > minfit):

			X.append(hag[i].Det1.time)
			X.append(hag[i].Det1.time)
#			X.append(hag[i].Det1.time)
			ypt=-1.*snr
			Y.append(ypt)
			Y.append(0.)
#			Y.append(snr)
			u.append(math.cos(hag[i].direction))
			v.append(math.sin(hag[i].direction))
			C.append('r')
			C.append('r')
			C.append('r')
			if(hag[i].direction_cc != NAN):		
				u.append(math.cos(hag[i].direction_cc))
				v.append(math.sin(hag[i].direction_cc))	
			else:
				u.append(0.)
				v.append(0.)		
			if(hag[i].direction_stk != NAN):
				u.append(math.cos(hag[i].direction_stk))
				v.append(math.sin(hag[i].direction_stk))
#			else:
#				u.append(0.)
#				v.append(0.)
			ang1 = hag[i].direction
			ang2 = hag[i].direction_cc
			ang3 = hag[i].direction_stk
			pee  = float(np.pi)

			ang1 = ang1-0.5*pee
			ang2 = ang2-0.5*pee
			ang3 = ang3-0.5*pee

			ang1 = ang1*(180./pee)
			ang2 = ang2*(180./pee)	
			ang3 = ang3*(180./pee)

#			print("Angles:", hag[i].direction, hag[i].direction_cc, ang1, ang2, pee)
			label = fm.format('{0}{1}',int(ang1),unichr(176))
			ax.text(hag[i].Det1.time,-0.8*snr,label)
			if(hag[i].direction_cc != NAN):	
				label = fm.format('{0}{1}',int(ang2),unichr(176))
				ax.text(hag[i].Det1.time,.2*snr,label)
			if(hag[i].direction_stk != NAN):
				label = fm.format('{0}{1}',int(ang3),unichr(176))
				ax.text(hag[i].Det1.time,1.2*snr,label)
			j += 1

#	plt.plot(pick1, snr2, 'ro')
#	plt.plot(pick1, snr3, 'mo')
#	plt.plot(pick1, snr1, 'bo')
	plt.plot(X, Y, 'mo')
	plt.quiver(X,Y,u,v,
#               C,
		units='height',
		scale=3.,
		width=0.03,
#		linewidth=(0.5,0.5),
		headaxislength=2., 
#		headlength=1.,
#		headwidth=2.,
		pivot='tail')
	label = fm.format('Direction of propagation')
	plt.xlabel(label)	
	
#	plt.ylabel("Direction of propagation")
	plt.ylim(-3*snr,3*snr)

	plt.subplot(514,sharex=ax1)
	power2=array('f',[])
	time1=array('f',[])
	for i in xrange(len(hag)):	
#		ppow2 = (math.log(hag[i].Det1.pow_30_50) + math.log(hag[i].Det2.pow_30_50) + math.log(hag[i].Det3.pow_30_50))/3.
		ppow2 = (hag[i].Det1.reclev + hag[i].Det2.reclev + hag[i].Det3.reclev)/3.
#		ppow2 = ppow2*20.
		power2.append(ppow2)
		time1.append(hag[i].Det1.time)
	
	if(len(power2) > 0):
		half = 0.5*(max(power2)+min(power2))
	
	for i in xrange(len(hag)):
		label = fm.format('{0}', i)
		plt.annotate(label,xy=(time1[i],power2[i]),xytext=(time1[i]-100., power2[i]-6.), 
			bbox=dict(boxstyle="square", fc="w"),arrowprops=dict(arrowstyle="->"))
		if(hag[i].whale_disc):
			plt.annotate('Whale',xy=(time1[i],power2[i]),xytext=(time1[i]-100., power2[i]+4.), 
			bbox=dict(boxstyle="round", fc="g"),arrowprops=dict(arrowstyle="->"))
		else:
			plt.annotate('N',xy=(time1[i],power2[i]),xytext=(time1[i]-100., power2[i]+4.), 
			bbox=dict(boxstyle="round", fc="r"),arrowprops=dict(arrowstyle="->"))

	label = fm.format('Log of amplitude and detection labelling')
	plt.xlabel(label)	
	plt.plot(time1,power2,'ro-')
	
	plt.subplot(515,sharex=ax1)
	cmap_gr=mp.colors.ListedColormap([[0.,0.,0.],[0.,0.,0.05],[0.,0.05,0.05],[0.05,0.05,0.05],[0.05,0.05,.1],[.05,0.1,.1],[.1,0.1,.1],[.1,0.1,.2],[.1,0.2,0.2],[0.2,0.2,0.2],[0.2,0.2,0.3],[.2,0.3,.3],[0.3,0.3,0.3],[.3,0.3,0.4],[1.,1.,0.],[1.,.0,.0],[1.,1.,1.]])
	Pxx,freqs,bins,im = pylab.specgram(wf1,NFFT=512,Fs=samprate,noverlap=511,cmap=cmap_gr)
#	plt.ylabel("Frequency in Hz.")

	label = fm.format('Seconds after {0}', pr_time)
	plt.xlabel(label)
#	plt.ylabel("Frequency in Hz.")
#	plt.subplots_adjust(bottom=0.2)

	plt.tight_layout(h_pad=.2)
	plt.ylim(1.,125.)
	plt.xlim(-10.,time_interval)

	pylab.savefig(figurename, bbox_inches='tight')
	plt.show()

