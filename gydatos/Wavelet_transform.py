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
['sta','Det_num','time','dt_cc','duration','freq','whale_disc','pow_15_25','pow_30_50','pow_50_60','spectim'], 
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
	

def MakeDetect(Det, w, snr, samprate, t_group, sta, shift):
	det = []
	spectim=[]
	init=0
	det_num=0
	for i in range(1,len(w)):
		t0=float(i)/samprate
		if(w[i] >= snr and w[i] > w[i-1] and init==0):
			init=1
			dw = (snr-w[i-1])/(w[i]-w[i-1])
			if dw != 0. :
				t = t0+shift+dw/samprate
			ndet = len(det)
#			print "ndet", ndet, det
			if((ndet==0) or (ndet > 0 and (t-det[ndet-1].time) > t_group)):
				newdet=[sta,det_num,t,0.,0.,0.,False,0.,0.,0.,spectim]
				det_num+=1
				print (t0, w[i-1], w[i])
				print (newdet,'\n')
				det.append(Detection._make(newdet))	
		if(w[i] < snr ):
			init=0	
	
	return det		



# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('wave_transform.cfg')

verbose = config.getint('control_parameters', 'verbose')
vector_plot = config.getboolean('control_parameters', 'vector_plot')

#wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

start_time  = config.get('time_interval','start_time')
end_time    = config.get('time_interval','end_time')
time_format = config.get('time_interval','format')
time_interval = config.getfloat('time_interval','time_interval') 
station1 = config.get('time_interval','Station') 

# station = config.get('time_interval','Station')

triad = config.get('Correlation','Triad')



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
start_level = config.getint('Filter','start_level')
wlt_type = config.get('Filter','wavelet_type')
wlt_side = config.get('Filter','wavelet_side')
add_levels = config.get('Filter','wavelet_add')
plot_wavelet = config.getboolean('Filter','plot_wavelet')

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

cc = array('f',[])
tt = array('f',[])
add = array('f',[])
dt1max, dt2max = gy.CompmaxDt(sfile,triad)

srate=gy.ReadSrate(wfdisc, station1, verbose)
print ("Sampling rate:", srate,'\n')
for i in range(1+int(time_interval*srate)):
		master_wf.append(0.)
samprate, start_master_wf, num_master = gy.ReadWfm(master_wf, wfdisc, start_master_e, end_master_e, sta_master, verbose)

for i in range(1+int(time_interval*samprate)):
		wf1.append(0.)
	
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
	
	#
	# Pad with zeroes until length reaches the next power of two
	#

	w = wlt.Wavelet('db1')
	N=1
	n=1
	while (len(wf1) > N):
		N = N*2
		if(len(wf1) < N):
			for i in range(N-len(wf1)):
				wf1.append(0.)
		len1 = wlt.swt_max_level(N)
		n1 = len(wf1)
		len3 = wlt.swt_max_level(n1)
		len2 = wlt.dwt_max_level(N, w.dec_len)
		print(N, 'swt:',len1,'dwt:',len2,len3,n1)
		n += 1
	for i in range(N):
		add.append(0.)	
	start_e += time_interval
	for i in range(N) :
		tt.append((float(i)/samprate))
	if (test_case):
		d1,d2 = gy.ReadSite(sfile, triad, 0)
		theta1 = (theta1+90.)*np.pi/180.
		dt1,dt2 = gy.CompDt(d1,d2,theta1,1.48)
		for i in xrange (len(wf1)):
			wf1[i] = noise_level*np.random.random_sample()
			
	#
	# Replace traces with boxcar in the middle of the trace 
	#
		for i in xrange(int(boxcar*samprate)):
			wf1[len(wf1)/2-int(boxcar*samprate/2)+i] += 1.
			
	
	if(whole_series_spectrum):
				
		label = fm.format('Wavelet Transform')
		print(label, len(wf1), wlt.swt_max_level(len(wf1)), start_level)
#		subf1.plot(freqs,20.*scipy.log10(fft1),'red')
#		[(cA2, cD2), (cA1, cD1)] = wlt.swt(wf1, 'db1', 2)
		[(cA7, cD7),(cA6, cD6),(cA5, cD5),(cA4, cD4),(cA3, cD3),(cA2, cD2), (cA1, cD1)] = wlt.swt(wf1, wlt_type, 7, start_level=start_level)
			
		if(plot_wavelet):
			figwvlt = plt.figure(figsize=(10,6))
			ax1=plt.subplot(711)
			plt.title(label)
			if(wlt_side == 'psi'):
				print('psi',start_level,int(add_levels[0]))
				plt.plot(tt,cD7)
				plt.subplot(712, sharex=ax1)
				plt.plot(tt,cD6)
				plt.subplot(713, sharex=ax1)
				plt.plot(tt,cD5)
				plt.subplot(714, sharex=ax1)
				plt.plot(tt,cD4)
				plt.subplot(715, sharex=ax1)
				plt.plot(tt,cD3)
				plt.subplot(716, sharex=ax1)
				plt.plot(tt,cD2)
				plt.subplot(717, sharex=ax1)
				plt.plot(tt,cD1)
			if(wlt_side == 'phi'):
				plt.plot(cA7)
				plt.subplot(712, sharex=ax1)
				plt.plot(cA6)
				plt.subplot(713, sharex=ax1)
				plt.plot(cA5)
				plt.subplot(714, sharex=ax1)
				plt.plot(cA4)
				plt.subplot(715, sharex=ax1)
				plt.plot(cA3)
				plt.subplot(716, sharex=ax1)
				plt.plot(cA2)
				plt.subplot(717, sharex=ax1)
				plt.plot(cA1)
			plt.show()
		if(len(add_levels) > 0):
			figadd = plt.figure(figsize=(10,6))
			ax1=plt.subplot(111)
			plt.title(label)
			for level in range(len(add_levels)):
				print("Level",add_levels[level])
				if(int(add_levels[level]) == 0):
					print("Adding cD7")		
					for i in range(N):
						add[i] += cD7[i]
				if(int(add_levels[level]) == 1):
					print("Adding cD6")		
					for i in range(N):
						add[i] += cD6[i]
				if(int(add_levels[level]) == 2):
					print("Adding cD5")		
					for i in range(N):
						add[i] += cD5[i]
				if(int(add_levels[level]) == 3):
					print("Adding cD4")		
					for i in range(N):
						add[i] += cD4[i]
				if(int(add_levels[level]) == 4):
					print("Adding cD3")		
					for i in range(N):
						add[i] += cD3[i]
				if(int(add_levels[level]) == 5):
					print("Adding cD2")		
					for i in range(N):
						add[i] += cD2[i]
				if(int(add_levels[level]) == 6):
					print("Adding cD1")		
					for i in range(N):
						add[i] += cD1[i]
			plt.plot(tt,add)
			plt.show()
	if(plot_crossco):
		if(envelope):
			hwf1 = scipy.signal.hilbert(cD1)
			print ("Hilbert 1")
#
# Envelope
#
			cD1 = np.absolute(hwf1)
			print ("Envelope")
		crossco0 = np.correlate(cD1,cD1,'full')
		print ("Correlate")
		crossco0 = np.append(crossco0, 0.)
		[(cA1, cD1)] = wlt.swt(crossco0, wlt_type, 1, start_level=4)
		trace = cA1
		trace_av = np.average(trace)
		sigi = 0.
		cross = 0.
		for i in range(len(trace)) :
				trace[i]=(trace[i]/trace_av)-1.
				if(time_int == 0):
					tt.append((float(i)/samprate)-len(trace)/samprate/2.)
#		for i in range(len(trace)/2) :
#				cross = cross + trace[i]*tt[i]
#				sigi = sigi + tt[i]*tt[i]
##		mom1 = cross/sigi
#		for i in range(len(trace)) :
#			trace[i] = trace[i]+2.*mom1*tt[i]
		i = 0

		start = 0
		
		if(time_int == 0):
			for i in range(len(trace)) :
				cc.append(time_int*0.5+trace[i])
		elif(time_int > 0):
			for i in range(len(trace)) :
				cc[i] = time_int*0.5+trace[i]
#		plt.plot([time_interval,(time_int-1)*0.5],[0.,time_int*0.5],'r',linewidth=5.)
		print(len(cc),len(tt))
		plt.plot(tt, cc, 'b')
#		plt.ylim(0.,10.)
		plt.xlim(0,time_interval/2.)
		print ("Crossco")
plt.show()
