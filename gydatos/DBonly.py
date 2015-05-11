#!/usr/bin/python

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
#import matplotlib.pyplot as plt
#import matplotlib.colors as cl
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
#import pylab
#import pyaudio
import wave
import glob
import os
import subprocess
import collections
import gydatos_ng as gy
import sqlite3 as sq

Detection = collections.namedtuple('Detection',
['sta','Det_num','time','dt_cc','duration','freq','whale_disc','pow_15_25','pow_30_50','pow_50_60','spectim'], 
verbose=True)
Hydro_Group = collections.namedtuple('HAG',['triad','HAG_num','whale_disc','Det1','Det2','Det3','direction','absfit','direction_cc','relfit','direction_stk','stkfit'], 
verbose=True)

fm = string.Formatter()


NAN=999.

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
# Look for a third detection in range and declare a group if fou
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
						whale = False
						newgrp = [triad,hag_num,whale,det1[i],det2[j],det3[k],OptDirection,maxdiff,NAN,0.,NAN,0.]
						
						print ("Dt vector ", dt)
						print ("New Group", maxdt1, maxdt2,i,j,k,OptDirection,'\n',newgrp,'\n')
						hag.append(Hydro_Group._make(newgrp))
						hag_num+=1
	return hag	
	
def GetId(con, object):
        get_string = "SELECT current from NEXTID where IDNAME='%s'"  %  (object)
        print(get_string)
	con.execute(get_string)
        data = int(con.fetchone()[0])
        print ("Latest id number for",object,"is: ",data)
        data=data+1
        get_string = "UPDATE NEXTID SET current=%d WHERE IDNAME='%s'"  %  (data, object)
        print(get_string)
        con.execute(get_string)
        return data
	
def MakeDetect(Det, w, snr, samprate, t_group, sta, shift, start_e):
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
				t = t0+shift+start_e+dw/samprate
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

time_interval  = config.getfloat('time_interval','time_interval')
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
#plot_detect = config.getboolean('Detect','plot_detect')
#plot_waveform = config.getboolean('Detect','plot_waveform')
#plot_crossco = config.getboolean('Detect','plot_crossco')
#stack = config.getboolean('Detect','stack')
#plot_stack = config.getboolean('Detect','plot_stack')
#plot_spect = config.getboolean('Detect','plot_spect')
#plot_hist = config.getboolean('Detect','plot_hist')
#plot_spectim = config.getboolean('Detect','plot_spectim')
#plot_power = config.getboolean('Detect','plot_power')
#plot_ratios = config.getboolean('Detect','plot_ratios')


sfile    = config.get('Azimuth','sfile')
halfwidth = config.getfloat('Azimuth','halfwidth')

refine_crossco = config.getboolean('Azimuth','refine_crossco')


if verbose > 2:
	print ('Order =', order, 'Low =', low, 'High =', high)

# Connect to database

con = sq.connect('hydro_whales.db')
cnx=con.cursor()
#cnx.execute("DELETE from WhaleDetection")
#cnx.execute("DELETE from  TRIAD")


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
	print ("Epoch end:", end_e,'\n')
	print ("Time interval:", time_interval)


dt1max, dt2max = gy.CompmaxDt(sfile,triad)

wf1=array('f',[])
wf2=array('f',[])
wf3=array('f',[])
master_wf=array('f',[])
srate=gy.ReadSrate(wfdisc, station1, verbose)
complete_interval = end_e - start_e
number_interval =int(complete_interval/time_interval)
for i in range(1+int(time_interval*srate)):
		wf1.append(0.)
for i in range(1+int(time_interval*srate)):
		wf2.append(0.)
for i in range(1+int(time_interval*srate)):
		wf3.append(0.)
	
Wn = np.array([])
bounds = np.array([float(low),float(high)])
Wn = np.append(Wn,bounds)
print(Wn)

b,a = scipy.signal.iirfilter(order, Wn)
print (a, b) 

StaLta = array('f',[])
for j in range(len(wf1)):
        StaLta.append(0.)

for interval in range(number_interval):

	start_int = start_e+interval*time_interval
	end_int = start_e+(interval+1)*time_interval

	samprate, start_wf1, num_samp1 = gy.ReadWfm(wf1, wfdisc, start_int, end_int, station1, verbose)
	samprate, start_wf2, num_samp2 = gy.ReadWfm(wf2, wfdisc, start_int, end_int, station2, verbose)
	samprate, start_wf3, num_samp3 = gy.ReadWfm(wf3, wfdisc, start_int, end_int, station3, verbose)
#
#
#
	print("Expected sample length:", int(1+time_interval*samprate))
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

	gy.CompStaLta(StaLta, L_ta, S_ta, samprate, fltwf1)
	det1=MakeDetect(Detection, StaLta, snr, samprate, t_group, station1, L_ta, start_int)

	gy.CompStaLta(StaLta, L_ta, S_ta, samprate, fltwf2)
	det2=MakeDetect(Detection, StaLta, snr, samprate, t_group, station2, L_ta, start_int)

	gy.CompStaLta(StaLta, L_ta, S_ta, samprate, fltwf3)
	det3=MakeDetect(Detection, StaLta, snr, samprate, t_group, station3, L_ta, start_int)

	hag=HydroGroupForm(sfile, triad, det1, det2, det3)
	gy.WhaleDiscriminant(hag,fltwf1,fltwf2,fltwf3,sfile,win_front,win_back,samprate,start_int)

	if(refine_crossco):
		gy.HydroGroupCCRefine(hag, fltwf1, fltwf2, fltwf3, sfile, triad, win_front, win_back,
		samprate, tcor, envelope, refine_crossco, halfwidth, start_int)
        
	del fltwf1
	del fltwf2
	del fltwf3
#
#
#

	Lshift = int(L_ta*samprate)
	Sshift = int(S_ta*samprate)
#
#
#
	for i in range(len(hag)):
                        detid1 = GetId(cnx,'detid')
                        insert_string = 'INSERT INTO WhaleDetection VALUES (%d,"%s",%f,%f,%d)' %  (detid1, hag[i].Det1.sta, hag[i].Det1.time,  hag[i].Det1.dt_cc, hag[i].Det1.whale_disc)
                        print insert_string
	 		cnx.execute(insert_string) 
                        detid2 = GetId(cnx,'detid')
                        insert_string = 'INSERT INTO WhaleDetection VALUES (%d,"%s",%f,%f,%d)' %  (detid2, hag[i].Det2.sta, hag[i].Det2.time,  hag[i].Det2.dt_cc, hag[i].Det2.whale_disc)
                        print insert_string
 	 		cnx.execute(insert_string)
                        detid3 = GetId(cnx,'detid')
                        insert_string = 'INSERT INTO WhaleDetection VALUES (%d,"%s",%f,%f,%d)' %  (detid3, hag[i].Det3.sta, hag[i].Det3.time,  hag[i].Det3.dt_cc, hag[i].Det3.whale_disc)
	 		cnx.execute(insert_string)
                        triadid = GetId(cnx,'triadid')
                        insert_string = 'INSERT INTO TRIAD VALUES ("%s",%d,%d,%d,%d,%f,%f,%f,%f,%f,%f)' %  (triad,triadid, detid1, detid2, detid3,  hag[i].direction, hag[i].absfit, hag[i].direction_cc, hag[i].relfit, hag[i].direction_stk, hag[i].stkfit)
                        print insert_string
	 		cnx.execute(insert_string) 

	con.commit()
if con:

	con.close()
