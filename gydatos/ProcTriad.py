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
import scipy.signal
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

Detection = collections.namedtuple('Detection',['sta','Det_num','time','dt_cc','duration','freq'], verbose=True)
Hydro_Group = collections.namedtuple('HAG',['triad','HAG_num','Det1','Det2','Det3','direction','direction_cc'], verbose=True)

fm = string.Formatter()

def HydroGroupCCRefine(hag,w1,w2,w3,win_front,win_back,samprate,tcor):
  '''
	Refine the delta-t times obtained from LTA/STA using crosscorrelation 
        and get an updated direction from these.
	'''	
#
# Compute either the crosscorrelation of the signal or the crosscorrelation of the envelopes of the signal in a window set to be win_front seconds before the STA-LTA pick and win_back seconds after.
#
	for i in range(len(hag)):		
#
# For each group, get the time window, envelope, and cross-correlation
#
		wf1 = w1[int((hag[i].Det1.time-win_front)*samprate):int((hag[i].Det1.time+win_back)*samprate)]
		wf2 = w2[int((hag[i].Det2.time-win_front)*samprate):int((hag[i].Det2.time+win_back)*samprate)]		
		wf3 = w3[int((hag[i].Det3.time-win_front)*samprate):int((hag[i].Det3.time+win_back)*samprate)]
#
# Hilbert transform
##
# Compute Cross-Correlations
#
		lcor = int(tcor*samprate)
		crossco0 = array('f',[])
		crossco1 = array('f',[])
		crossco2 = array('f',[])
		crossco3 = array('f',[])
		crossco0 = np.correlate(wf1,wf1,'full')
		crossco1 = np.correlate(wf1,wf2,'full')
		#
		#
		#
		print ("Crosscorrelation 1")
		crossco2 = np.correlate(wf1,wf3,'full')
		print ("Crosscorrelation 2")

		if(envelope):
			hwf1 = scipy.signal.hilbert(crossco0)
			print ("Hilbert 1")
			hwf2 = scipy.signal.hilbert(crossco1)
			print ("Hilbert 2")
			hwf3 = scipy.signal.hilbert(crossco2)
			print ("Hilbert 3")
#
# Envelope
#
			wf1 = np.absolute(hwf1)
			wf2 = np.absolute(hwf2)
			wf3 = np.absolute(hwf3)

#
# Compute Cross-Correlations
#
		t0 = gy.TimeOfMax(wf1)
		t1 = gy.TimeOfMax(wf2)
		t2 = gy.TimeOfMax(wf3)
		dt1 = (t1-t0)/samprate
		dt2 = (t2-t0)/samprate
		print ("HAG number, HAG:", i, hag[i])
		
		Detnew = hag[i].Det1._replace(dt_cc=hag[i].Det1.time)
		hag[i] = hag[i]._replace(Det1=Detnew)
		Detnew = hag[i].Det2._replace(dt_cc=hag[i].Det2.time+dt1)
		hag[i] = hag[i]._replace(Det2=Detnew)		
		Detnew = hag[i].Det3._replace(dt_cc=hag[i].Det3.time+dt2)
		hag[i] = hag[i]._replace(Det3=Detnew)

		dt = (hag[i].Det2.dt_cc - hag[i].Det1.dt_cc), (hag[i].Det3.dt_cc - hag[i].Det1.dt_cc)
						
		OptDirection, theta, radii = gy.CompOptDirection(dt,sfile,triad)
		print ("New direction:", OptDirection)
		hag[i]=hag[i]._replace(direction_cc=OptDirection)
		print ("HAG number, New HAG:", i, hag[i])

def HydroGroupForm(sfile, triad, det1, det2, det3):
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
						dt=dt1,dt2
						
						OptDirection,theta,radii = gy.CompOptDirection(dt,sfile,triad)
						fig = plt.figure(figsize=(6,6))
						ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
						N = len(radii)
						width = 2.*np.pi/N
						bars = ax.bar(theta, radii, width=width, bottom=0.)
						for r,bar in zip(radii, bars):
							bar.set_facecolor( cm.jet(r/10.))
							bar.set_alpha(0.5)
	
						newgrp = [triad,hag_num,det1[i],det2[j],det3[k],OptDirection,0.]
							
						label = fm.format('HAG nb {0}', hag_num)
						label = triad+' '+label

						plt.title(label)
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
				newdet=[sta,det_num,t,0.,0.,0.]
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

sfile    = config.get('Azimuth','sfile')
az_width = config.getfloat('Azimuth','width')
az_type  = config.get('Azimuth','type')
test_case = config.getboolean('Azimuth','test_case')


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


dt1max, dt2max = gy.CompmaxDt(sfile,triad)
	

samprate, start_wf1, wf1 = gy.ReadWfm(wfdisc, start_e, end_e, station1, verbose)
samprate, start_wf2, wf2 = gy.ReadWfm(wfdisc, start_e, end_e, station2, verbose)
samprate, start_wf3, wf3 = gy.ReadWfm(wfdisc, start_e, end_e, station3, verbose)
samprate, start_master_wf, master_wf = gy.ReadWfm(wfdisc, start_master_e, end_master_e, sta_master, verbose)
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
if envelope:
	hwf1 = scipy.signal.hilbert(fltwf1)
	print ("Hilbert 1")
	hwf2 = scipy.signal.hilbert(fltwf2)
	print ("Hilbert 2")
	hwf3 = scipy.signal.hilbert(fltwf3)
	print ("Hilbert 3")
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
if (correlate):
	lcor=int(tcor*samprate)
	crossco0 = array('f',[])
	crossco1 = array('f',[])
	crossco2 = array('f',[])
	crossco3 = array('f',[])
	crossco0 = np.correlate(fltwf1,fltwf1,'full')
	crossco1 = np.correlate(fltwf1,fltwf2,'full')
	#
	#
	#
	print ("Crosscorrelation 1")
	crossco2 = np.correlate(fltwf1,fltwf3,'full')
	print ("Crosscorrelation 2")
	t0 = gy.TimeOfMax(crossco0)
	t1 = gy.TimeOfMax(crossco1)
	t2 = gy.TimeOfMax(crossco2)
	dt1 = (t1-t0)/samprate
	dt2 = (t2-t0)/samprate
	crossco3 = np.correlate(fltwf2,fltwf3,'full')
	print ("Crosscorrelation 3")
	t3 = gy.TimeOfMax(crossco3)
	dt3= (t3-t0)/samprate
	diffdt = dt3-(dt2+dt1)
	print ("t0, t1, t2, t3", t0, t1, t2, t3)
	print ("dt1, dt2, dt3, diffdt", dt1, dt2, dt3, diffdt)

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

hag=HydroGroupForm(sfile, triad, det1, det2, det3)
HydroGroupCCRefine(hag,fltwf1,fltwf2,fltwf3,win_front,win_back,samprate,tcor)

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
	for i in range(len(hag)):
		pick1.append(float(Lshift-Sshift)/samprate+hag[i].Det1.time)		
		pick2.append(float(Lshift-Sshift)/samprate+hag[i].Det2.time)		
		pick3.append(float(Lshift-Sshift)/samprate+hag[i].Det3.time)
		snr1.append(snr)
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
	
	for i in range(len(hag)):		
		snr2.append(-1.*snr)

	for i in range(len(det1)):
		allpick1.append(det1[i].time)
		tpick1.append(0.)
	
	plt.plot(tflt,fltwf1, 'black',allpick1,tpick1, 'rs', pick1, snr2, 'ro')
		
	label = fm.format(' {0} - Filtered',station1)
	plt.ylabel(label)
	ax=plt.subplot(313,sharex=ax1)
        
	X = array('f',[])
	Y = array('f',[])
	u = array('f',[])
	v = array('f',[])
	for i in range(len(hag)):
		X.append(pick1[i])
		X.append(pick1[i])
		Y.append(-1.*snr)
		Y.append(snr)
		u.append(math.cos(hag[i].direction))
		v.append(math.sin(hag[i].direction))		
		u.append(math.cos(hag[i].direction_cc))
		v.append(math.sin(hag[i].direction_cc))
		label = fm.format('Angle {0} ',180.*(hag[i].direction-0.5*np.pi)/np.pi)
		ax.text(pick1[i],-0.8*snr,label)
		label = fm.format('Angle {0} ',180.*(hag[i].direction_cc-0.5*np.pi)/np.pi)
		ax.text(pick1[i],1.2*snr,label)

	plt.plot(pick1, snr2, 'ro')
	plt.plot(pick1, snr1, 'bo')
	plt.quiver(X,Y,u,v)
	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)	
	
	plt.ylabel("Incoming direction of signal")
	plt.ylim(-2*snr,2*snr)
	plt.show()
else:
	t=array('f',[])
	for j in range (len(fltwf1)):
		t.append(float(j)/samprate)		
	print ("samprate = ", samprate, "length =", len(fltwf1))
	print ("Median1=", median1, "Median2=", median2, "Median3=", median3)
	ax1=plt.subplot(331)
	plt.plot(t, fltwf1, 'b', tmed, med1,'r')
	
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



	pick1 = array('f',[])
	snr1 = array('f',[])
	for i in range(len(det1)):
		pick1.append(float(Lshift-Sshift)/samprate+det1[i].time)
		snr1.append(snr)
	plt.plot(tsta,StaLta1,pick1,snr1, 'o')
	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)
	plt.ylabel("StaLta ratio")

	plt.subplot(335,sharex=ax1)
	pick2 = array('f',[])
	snr2 = array('f',[])
	for i in range(len(det2)):
		pick2.append(float(Lshift-Sshift)/samprate+det2[i].time)
		snr2.append(snr)
	
	plt.plot(tsta,StaLta2, pick2, snr2, 'o')
	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)
	plt.ylabel("StaLta ratio")

	plt.subplot(338,sharex=ax1)

	pick3 = array('f',[])
	snr3 = array('f',[])
	for i in range(len(det3)):
		pick3.append(float(Lshift-Sshift)/samprate+det3[i].time)
		snr3.append(snr)

	plt.plot(tsta, StaLta3, pick3, snr3, 'o')
	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)
	plt.ylabel("StaLta ratio")
	plt.subplot(333,sharex=ax1)

	
	if (correlate):
		tcor=array('f',[])
		for j in range (len(crossco1)):
			tcor.append(-0.5*float(len(crossco1)-1)/samprate+float(j)/samprate)	
		plt.plot(tcor,crossco1)
	elif (plot_detect):
		plt.plot(pick1,snr1,'go',pick2,snr2,'rx',pick3,snr3,'bs')
	else:
		plt.plot(tsta,StaLta1)

	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)
	plt.ylabel("StaLta ratio")
	plt.subplot(336,sharex=ax1)

	if (correlate):	
		plt.plot(tcor,crossco2)
	else:
		plt.plot(tsta,StaLta2)	

	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)
	plt.ylabel("StaLta ratio")
	plt.subplot(339,sharex=ax1)
	
	if (correlate):	
		plt.plot(tcor,crossco3)
	else:
		plt.plot(tsta,StaLta3)

	label = fm.format('Seconds after {0}',start_time)
	plt.xlabel(label)
	plt.ylabel("StaLta ratio")

	#
	#plt.subplots_adjust(bottom=0.2)
	#
	figurename = triad+'_'+start_time+'.png'
	pylab.savefig(figurename)
	plt.show()
