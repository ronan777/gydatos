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
import matplotlib.mpl as mp
import matplotlib.colorbar as col
import matplotlib.pyplot as plt
import matplotlib.colors as cl
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
import _mysql as ms

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
config.read('Waveforms.cfg')

verbose = config.getint('control_parameters', 'verbose')
vector_plot = config.getboolean('control_parameters', 'vector_plot')

#wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

start_time  = config.get('time_interval','start_time')
end_time    = config.get('time_interval','end_time')
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

# Connect to database

cnx = ms.connect(host='sql4.freesqldatabase.com',user='sql457121', passwd='fN2%aL4!', db='sql457121')

cnx.query('TRUNCATE table WhaleDetection')

print ('Connector:', cnx)
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

wf1=array('f',[])
wf2=array('f',[])
wf3=array('f',[])
master_wf=array('f',[])
srate=gy.ReadSrate(wfdisc, station1, verbose)
time_interval = end_e - start_e
for i in range(1+int(time_interval*srate)):
		wf1.append(0.)
for i in range(1+int(time_interval*srate)):
		wf2.append(0.)
for i in range(1+int(time_interval*srate)):
		wf3.append(0.)
for i in range(1+int(time_interval*srate)):
		master_wf.append(0.)
	

samprate, start_wf1, num_samp1 = gy.ReadWfm(wf1, wfdisc, start_e, end_e, station1, verbose)
samprate, start_wf2, num_samp2 = gy.ReadWfm(wf2, wfdisc, start_e, end_e, station2, verbose)
samprate, start_wf3, num_samp3 = gy.ReadWfm(wf3, wfdisc, start_e, end_e, station3, verbose)
samprate, start_master_wf, num_master = gy.ReadWfm(master_wf, wfdisc, start_master_e, end_master_e, sta_master, verbose)
#
#
#
if (test_case):
	d1,d2 = gy.ReadSite(sfile, triad, 0)
	theta1 = (theta1+90.)*np.pi/180.
	dt1,dt2 = gy.CompDt(d1,d2,theta1,1.48)
	for i in range (len(wf1)):
		wf1[i] = noise_level*np.random.random_sample()
		wf2[i] = noise_level*np.random.random_sample()
		wf3[i] = noise_level*np.random.random_sample()
#
# Replace traces with boxcar in the middle of the trace 
#

	for i in range(int(boxcar*samprate)):
		wf1[len(wf1)/2-int(boxcar*samprate/2)+i] += 1.
		wf2[len(wf2)/2-int(boxcar*samprate/2)+i+int(dt1*samprate)] += 1.
		wf3[len(wf3)/2-int(boxcar*samprate/2)+i+int(dt2*samprate)] += 1.


Wn = np.array([])
bounds = np.array([float(low),float(high)])
Wn = np.append(Wn,bounds)
print(Wn)

b,a = scipy.signal.iirfilter(order, Wn)
print (a, b) 
print("Expected sample length:", int(1+(end_e-start_e)*samprate))
print("Read", num_samp1, num_samp2, num_samp3, "samples")
fltwf1 = scipy.signal.lfilter(b, a, wf1)
print ("Filtered 1, length:", len(wf1))
fltwf2 = scipy.signal.lfilter(b, a, wf2)
print ("Filtered 2, length:", len(wf2))
fltwf3 = scipy.signal.lfilter(b, a, wf3)
print ("Filtered 3, length:", len(wf3))

#
# STA/LTA detector 
#
#
# Detections based on STA/LTA 
#	

StaLta1=gy.CompStaLta(L_ta, S_ta, samprate, fltwf1)
det1=MakeDetect(Detection, StaLta1, snr, samprate, t_group, station1, L_ta)

StaLta2=gy.CompStaLta(L_ta, S_ta, samprate, fltwf2)
det2=MakeDetect(Detection, StaLta2, snr, samprate, t_group, station2, L_ta)

StaLta3=gy.CompStaLta(L_ta, S_ta, samprate, fltwf3)
det3=MakeDetect(Detection, StaLta3, snr, samprate, t_group, station3, L_ta)


hag=HydroGroupForm(sfile, triad, det1, det2, det3, plot_detect)
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
# Plotting of filtered waveforms
#

if(vector_plot == True):	

	tsta  = np.array([])
	tsta2 = np.array([])
	tsta3 = np.array([])

	for j in range (len(StaLta1)):
		tsta = np.append(tsta,[float(Lshift)/samprate+float(j)/samprate])	
	
	ax1=plt.subplot(512)

	pick1 = np.array([])	
	pick2 = np.array([])	
	pick3 = np.array([])
	snr1 = np.array([])
	detid=0
	nwhale = 0
	for i in range(len(hag)):
#		if(hag[i].whale_disc):
			pick1 = np.append(pick1,[hag[i].Det1.time])		
			pick2 = np.append(pick2,[hag[i].Det2.time])		
			pick3 = np.append(pick3,[hag[i].Det3.time])
#	 		cnx.query("""INSERT INTO WhaleDetection VALUES (detid, hag[i].Det1.sta, hag[i].Det1.time, hag[i].Det1.whale_disc)""")
                        insert_string = 'INSERT INTO WhaleDetection VALUES (%d,"%s",%f,%f,%d)' %  (detid, hag[i].Det1.sta, hag[i].Det1.time,  hag[i].Det1.dt_cc, hag[i].Det1.whale_disc)
	 		cnx.query(insert_string)                      
                        insert_string = 'INSERT INTO WhaleDetection VALUES (%d,"%s",%f,%f,%d)' %  (detid, hag[i].Det2.sta, hag[i].Det2.time,  hag[i].Det2.dt_cc, hag[i].Det2.whale_disc)
 	 		cnx.query(insert_string)
                        insert_string = 'INSERT INTO WhaleDetection VALUES (%d,"%s",%f,%f,%d)' %  (detid, hag[i].Det3.sta, hag[i].Det3.time,  hag[i].Det3.dt_cc, hag[i].Det3.whale_disc)
	 		cnx.query(insert_string)
			snr1 = np.append(snr1,[snr])
			nwhale += 1
	print(len(tsta),len(StaLta1),len(pick1),len(snr1))
        cnx.close()
	plt.plot(tsta,StaLta1,'red',  pick1,snr1, 'ro')

#	print(len(tsta),len(StaLta1),len(StaLta2),len(StaLta3))	

	plt.plot(tsta,StaLta2,'green',pick2,snr1, 'go')
	plt.plot(tsta,StaLta3,'blue' ,pick3,snr1, 'bo')
	plt.show()
	plt.xlabel("STA/LTA ratio")

	plt.subplot(511,sharex=ax1)
	
	snr2 = np.array([])
	snr3 = np.array([])
	snr1 = np.array([])
	tflt = np.array([])
	tpick1 = np.array([])
	allpick1 = np.array([])
	for j in range (len(fltwf1)):
		tflt = np.append(tflt,[float(j)/samprate])		
	hmax = 0.5*max(fltwf1)
	
	for i in range(len(pick1)):
		snr1=np.append(snr1,[-1.*hmax])
	for i in range(len(pick2)):
		snr2=np.append(snr2,[0.])
	for i in range(len(pick3)):
		snr3=np.append(snr3,[0.])
		
	
	plt.plot(tflt,fltwf1, 'black',pick1,snr1, 'ro', pick2, snr2, 'ys', pick3, snr3, 'ys')
		
	label = fm.format(' {0} - Filtered',station1)
	plt.xlabel(label)
	ax=plt.subplot(513,sharex=ax1)
        
	X = np.array([])
	Y = np.array([])
	u = np.array([])
	v = np.array([])
	j = 0

	for i in range(len(hag)):
		if(hag[i].whale_disc and (hag[i].absfit > minfit or hag[i].relfit > minfit)):
#		if(hag[i].absfit > minfit or hag[i].relfit > minfit):

			X = np.append(X,[hag[i].Det1.time])
			X = np.append(X,[hag[i].Det1.time])
			X = np.append(X,[hag[i].Det1.time])
			ypt=-1.*snr
			Y = np.append(Y,[ypt])
			Y = np.append(Y,[0.])
			Y = np.append(Y,[snr])
			u = np.append(u,[math.cos(hag[i].direction)])
			v = np.append(v,[math.sin(hag[i].direction)])
			if(hag[i].direction_cc != NAN):		
				u = np.append(u, [math.cos(hag[i].direction_cc)])
				v = np.append(v, [math.sin(hag[i].direction_cc)])	
			else:
				u = np.append(u, [0.])
				v = np.append(v, [0.])		
			if(hag[i].direction_stk != NAN):
				u = np.append(u, [math.cos(hag[i].direction_stk)])
				v = np.append(v, [math.sin(hag[i].direction_stk)])
			else:
				u = np.append(u, [0.])
				v = np.append(v, [0.])
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

	plt.plot(pick1, snr2, 'ro')
	plt.plot(pick1, snr3, 'mo')
	plt.plot(pick1, snr1, 'bo')
	plt.quiver(X,Y,u,v,units='height',scale=1.5)
	label = fm.format('Direction of propagation')
	plt.xlabel(label)	
	
#	plt.ylabel("Direction of propagation")
	plt.ylim(-3*snr,3*snr)

	plt.subplot(514,sharex=ax1)
	power2=np.array([])
	time1=np.array([])
	for i in range(len(hag)):	
		ppow2 = (math.log(hag[i].Det1.pow_30_50) + math.log(hag[i].Det2.pow_30_50) + math.log(hag[i].Det3.pow_30_50))/3.
		power2 = np.append(power2, [ppow2])
		time1 = np.append(time1,[hag[i].Det1.time])
	if(len(power2) == 0):
		power2 = np.append(power2,[0.])
		time1 = np.append(time1,[0.])
	
	half = 0.5*(max(power2)+min(power2))
	
	for i in range(len(hag)):
		label = fm.format('{0}', i)
		plt.annotate(label,xy=(time1[i],power2[i]),xytext=(time1[i]-100., power2[i]-.3), 
			bbox=dict(boxstyle="square", fc="w"),arrowprops=dict(arrowstyle="->"))
		if(hag[i].whale_disc):
			plt.annotate('Whale',xy=(time1[i],power2[i]),xytext=(time1[i]-100., power2[i]+.2), 
			bbox=dict(boxstyle="round", fc="g"),arrowprops=dict(arrowstyle="->"))
		else:
			plt.annotate('N',xy=(time1[i],power2[i]),xytext=(time1[i]-100., power2[i]+.2), 
			bbox=dict(boxstyle="round", fc="r"),arrowprops=dict(arrowstyle="->"))

	plt.plot(time1,power2,'ro-')	
	
	cmap_gr = mp.colors.ListedColormap([[0.,0.,0.],[0.,0.,0.05],[0.,0.05,0.05],[0.05,0.05,0.05],[0.05,0.05,.1],[.05,0.1,.1],[.1,0.1,.1],
                                          [.1,0.1,.2],[.1,0.2,0.2],[0.2,0.2,0.2],[0.2,0.2,0.3],[.2,0.3,.3],[0.3,0.3,0.3],[.3,0.3,0.4],[1.,1.,.0],[1.,.0,.0],[1.,1.,1.]])
	Pxx,freqs,bins,im = pylab.specgram(wf1,NFFT=512,Fs=samprate,noverlap=511,cmap=cmap_gr)
	print('Pxx.shape',Pxx.shape[0],freqs.shape[0],bins.shape[0])
        cnx = ms.connect(host='sql4.freesqldatabase.com',user='sql457121', passwd='fN2%aL4!', db='sql457121')
        cnx.query("SELECT * from WhaleDetection")
	result=cnx.store_result()
        
	print result.fetch_row(6)
        cnx.close()

#	plt.ylabel("Frequency in Hz.")
	sp=plt.subplot(515,sharex=ax1)
	plt.pcolor( Pxx)
	label = fm.format('Seconds after {0}', start_time)
	plt.xlabel(label)
#	plt.ylabel("Frequency in Hz.")
#	plt.subplots_adjust(bottom=0.2)

	plt.tight_layout(h_pad=.2)
	plt.ylim(15.,50.)
	plt.xlim(-60.,tflt[len(tflt)-1]+60.)

#	cl.LogNorm(vmin=0., vmax=50., clip=True)
	cb=plt.colorbar()
	cb.ax.set_yticks([-100., 100])
	plt.show()

if(plot_hist):
	theta=np.linspace(0.,2.*np.pi,37)

	hist_abs=array('f',[])
	hist_rel=array('f',[])
	hist_stk=array('f',[])
	for i in range(37):
		hist_abs.append(0.)
		hist_rel.append(0.)
		hist_stk.append(0.)
	
	nwhale1 = 0
	nwhale2 = 0
	nwhale3 = 0
	for i in range(len(hag)):
		if(hag[i].whale_disc): 
			if (hag[i].absfit >= minfit):			
				hist_abs[int(18.*(hag[i].direction/np.pi))] += 1.
				nwhale1 += 1
			if (hag[i].relfit >= minfit):
				hist_rel[int(18.*(hag[i].direction_cc/np.pi))] += 1.	
				nwhale2 += 1
			if (hag[i].stkfit >= minfit and hag[i].direction_stk != NAN):
				hist_stk[int(18.*(hag[i].direction_stk/np.pi))] += 1.
				nwhale3 += 1

	fig1 = plt.figure(figsize=(6,6))
	ax = fig1.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
	N = len(theta)
	width = 2.*np.pi/N
	bars = ax.bar(theta, hist_abs, width=width, bottom=0.)
	for r,bar in zip(hist_abs, bars):
		bar.set_facecolor( cm.jet(r/10.))
		bar.set_alpha(0.6)
		
	label = fm.format('Histogram abs {0} points', nwhale1)
	label = triad+' '+label
	plt.title(label)

	fig2 = plt.figure(figsize=(6,6))
	ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
	N = len(theta)
	width = 2.*np.pi/N
	bars = ax2.bar(theta, hist_rel, width=width, bottom=0.)
	for r,bar in zip(hist_rel, bars):
		bar.set_facecolor( cm.jet(r/10.))
		bar.set_alpha(0.6)
		
	label = fm.format('Histogram rel {0} points', nwhale2)
	label = triad+' '+label
	plt.title(label)

	fig3 = plt.figure(figsize=(6,6))
	ax3 = fig3.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
	N = len(theta)
	width = 2.*np.pi/N
	bars = ax3.bar(theta, hist_stk, width=width, bottom=0.)
	for r,bar in zip(hist_stk, bars):
		bar.set_facecolor( cm.jet(r/10.))
		bar.set_alpha(0.6)
		
	label = fm.format('Histogram stk {0} points', nwhale3)
	label = triad+' '+label
	plt.title(label)
	plt.show()

if(plot_power):
	power1=array('f',[])
	power2=array('f',[])
	power3=array('f',[])
	time1=array('f',[])
	time2=array('f',[])
	time3=array('f',[])
	for i in range(len(hag)):
		ppow1 = (math.log(hag[i].Det1.pow_15_25) + math.log(hag[i].Det2.pow_15_25) + math.log(hag[i].Det3.pow_15_25))/3.
		ppow2 = (math.log(hag[i].Det1.pow_30_50) + math.log(hag[i].Det2.pow_30_50) + math.log(hag[i].Det3.pow_30_50))/3.
		ppow3 = (math.log(hag[i].Det1.pow_50_60) + math.log(hag[i].Det2.pow_50_60) + math.log(hag[i].Det3.pow_50_60))/3.

		power1.append(ppow1)
		power2.append(ppow2)
		power3.append(ppow3)
		time1.append(hag[i].Det1.time)
		
	plt.plot(time1,power1, 'ro',time1,power2,'go',time1,power3,'bo')
	plt.show()

if(plot_ratios):
#	power1=array('f',[])
	power2=array('f',[])
	power3=array('f',[])
	time1=array('f',[])
#	time2=array('f',[])
#	time3=array('f',[])
	for i in range(len(hag)):
		ppow1 = (math.log(hag[i].Det1.pow_15_25) + math.log(hag[i].Det2.pow_15_25) + math.log(hag[i].Det3.pow_15_25))/3.
		ppow2 = (math.log(hag[i].Det1.pow_30_50) + math.log(hag[i].Det2.pow_30_50) + math.log(hag[i].Det3.pow_30_50))/3.
		ppow3 = (math.log(hag[i].Det1.pow_50_60) + math.log(hag[i].Det2.pow_50_60) + math.log(hag[i].Det3.pow_50_60))/3.

#		power1.append(ppow1)
		power2.append(ppow2/ppow1)
		power3.append(ppow2/ppow3)
		time1.append(hag[i].Det1.time)
		
	plt.xlim(-60.,time1[len(time1)-1]+60.)
	plt.plot(time1,power2,'go',time1,power3,'bo')
	plt.show()	

