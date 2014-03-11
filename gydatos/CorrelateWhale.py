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
import matplotlib.patches as mpatches 
from   matplotlib.widgets import Button, Slider
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
import gydatos as gy

epsilon = 0.0001



# Read configuration parameters

config = ConfigParser.RawConfigParser()
config.read('Waveforms.cfg')

verbose = config.getint('control_parameters', 'verbose')

#wfile = config.get('waveforms', 'wfile')
wfdisc = config.get('waveforms', 'wfdisc')

envelope = config.getboolean('Filter','envelope')
order = config.getint('Filter','order')
low = config.get('Filter','low')
high = config.get('Filter','high')

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
halfwidth = config.get('Azimuth','halfwidth')

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

samprate, null, wf1 = gy.ReadWfm(wfdisc, start_e, end_e, station1, verbose)
samprate, null, wf2 = gy.ReadWfm(wfdisc, start_e, end_e, station2, verbose)
samprate, null, wf3 = gy.ReadWfm(wfdisc, start_e, end_e, station3, verbose)
samprate, null, wfm = gy.ReadWfm(wfdisc, start_master_e, end_master_e, sta_master, verbose)


# Design digital IIR filter

Wn = array('f',[])
Wn.append(float(low))
Wn.append(float(high))
b,a = scipy.signal.iirfilter(order, Wn)
print (a, b) 

wf1 = scipy.signal.lfilter(b, a, wf1)
wf2 = scipy.signal.lfilter(b, a, wf2)
wf3 = scipy.signal.lfilter(b, a, wf3)
wf4 = scipy.signal.lfilter(b, a, wf3)
wfm = scipy.signal.lfilter(b, a, wfm)

# Calculate cross-correlations of master with all three waveforms


cross1 = np.correlate(wf1, wfm,'full')
cross2 = np.correlate(wf2, wfm,'full')
cross3 = np.correlate(wf3, wfm,'full')
cross4 = np.correlate(wf2, wf3,'full')



#
# Plotting of autocorrelations
#
m0=cross1.max()
t=array('f',[])
for j in range (len(cross1)):
	t.append(float(j)/samprate)
	cross1[j] = cross1[j]/m0
	cross2[j] = cross2[j]/m0
	cross3[j] = cross3[j]/m0
	cross4[j] = cross4[j]/m0

		
if envelope:
	hwf1 = scipy.signal.hilbert(cross1)
	print "Hilbert 1"
	hwf2 = scipy.signal.hilbert(cross2)
	print "Hilbert 2"
	hwf3 = scipy.signal.hilbert(cross3)
	print "Hilbert 3"
	hwf4 = scipy.signal.hilbert(cross4)
	wf1 = np.absolute(hwf1)
	wf2 = np.absolute(hwf2)
	wf3 = np.absolute(hwf3)	
	wf4 = np.absolute(hwf4)

twf=array('f',[])
maxtime=array('f',[])
maxvalue=array('f',[])
for j in range (len(wf1)):
	twf.append(float(j)/samprate)	
print "samprate = ", samprate, "length =", len(cross1), "Maximum of CC:", m0
maxtime.append(120.)
maxvalue.append(2.)

ax1=plt.subplot(421)
ax1.set_xlim([0.,0.5])
plt.xlim(70.,130.)
w1, = plt.plot(t, cross1)
fm = string.Formatter()
label = fm.format(' {0} - {1} - Crosscorrelation',station1,station1)
plt.title(label)

ax1=plt.subplot(422, sharex=ax1)
maxvalue[0] = wf1.max()
sumvalue = sum(wf1)/samprate
for j in range (len(wf1)):
	if(maxvalue[0]-epsilon <= wf1[j]):
		maxtime[0] = float(j);
for j in range (len(wf1)):
	wf1[j] = wf1[j]/sumvalue

maxtime[0] = maxtime[0]/samprate
maxvalue[0] = maxvalue[0]/sumvalue

#
# Fit a Laplace distribution or a Gaussian
#
halfwidth = float(halfwidth)
lapx0, lapstdv  ,sumd = gy.FitLaplace(wf1,samprate, halfwidth)
lapdist = array('f',[])
gaudist = array('f',[])
for j in range(len(twf)):
#	print maxvalue[0], lapx0, lapstdv
	arg = 0.5*(1./lapstdv)*math.exp(-1.*abs(j/samprate-lapx0)/lapstdv)
#	print(arg)
	lapdist.append(sumd*arg)
	arg = (1./math.sqrt(2.*math.pi))*(1./lapstdv)*math.exp(-1.*(j/samprate-lapx0)*(j/samprate-lapx0)/2./lapstdv/lapstdv)
	gaudist.append(sumd*arg)

maxtime.append(lapx0)
maxvalue.append(max(lapdist))

plt.plot(twf, wf1,'blue', maxtime, maxvalue,'ro', twf, lapdist, 'red', twf, gaudist, 'green')

label = fm.format(' {0} - {1} - CC Envelope',station1,station1)
plt.title(label)

plt.subplot(423, sharex=ax1)
w2, = plt.plot(t, cross2)
label = fm.format(' {0} - {1} - Crosscorrelation',station1,station2)
plt.title(label)


plt.subplot(424, sharex=ax1)
maxvalue[0] = wf2.max()
sumvalue=sum(wf2)/samprate
for j in range (len(wf2)):
	if(maxvalue[0]-epsilon <= wf2[j]):
		maxtime[0] = float(j);
for j in range (len(wf2)):
	wf2[j] = wf2[j]/sumvalue		

maxtime[0] = maxtime[0]/samprate
maxvalue[0] = maxvalue[0]/sumvalue

lapx0, lapstdv, sumd = gy.FitLaplace(wf2,samprate,halfwidth)
maxtime[1] = lapx0


for j in range(len(twf)):
#	print maxvalue[0], lapx0, lapstdv
	arg = 0.5*(1./lapstdv)*math.exp(-1.*abs(j/samprate-lapx0)/lapstdv)
#	print(arg)
	lapdist[j]=sumd*arg
	arg = (1./math.sqrt(2.*math.pi))*(1./lapstdv)*math.exp(-1.*(j/samprate-lapx0)*(j/samprate-lapx0)/2./lapstdv/lapstdv)
	gaudist[j]=sumd*arg

maxvalue[1] = max(lapdist)

plt.plot(twf, wf2, 'blue', maxtime, maxvalue,'ro', twf, lapdist, 'red', twf, gaudist, 'green')

label = fm.format(' {0} - {1} - CC Envelope',station1, station2)
plt.title(label)

ax1=plt.subplot(425, sharex=ax1)
w3, = plt.plot(t, cross3)
label = fm.format(' {0} - {1} - Crosscorrelation',station1, station3)
plt.title(label)

plt.ylabel("Correlation coefficient")
plt.title(label)
plt.subplot(426, sharex=ax1)
maxvalue[0] = wf3.max()
sumvalue = sum(wf3)/samprate
for j in range (len(wf3)):
	if(maxvalue[0]-epsilon <= wf3[j]):
		maxtime[0] = float(j);
for j in range (len(wf3)):
	wf3[j] = wf3[j]/sumvalue
maxtime[0] = maxtime[0]/samprate
maxvalue[0] = maxvalue[0]/sumvalue


lapx0, lapstdv , sumd = gy.FitLaplace(wf3,samprate, halfwidth)
maxtime[1] = lapx0


for j in range(len(twf)):
	arg = 0.5*(1./lapstdv)*math.exp(-1.*abs(j/samprate-lapx0)/lapstdv)
	lapdist[j] = sumd*arg
	arg = (1./math.sqrt(2.*math.pi))*(1./lapstdv)*math.exp(-1.*(j/samprate-lapx0)*(j/samprate-lapx0)/2./lapstdv/lapstdv)
	gaudist[j] = sumd*arg
maxvalue[1] = max(lapdist)
plt.plot(twf, wf3, 'blue', maxtime, maxvalue,'ro', twf, lapdist, 'red', twf, gaudist, 'green')

label = fm.format(' {0} - {1} - CC Envelope',station1,station3)
plt.title(label)

ax1=plt.subplot(427, sharex=ax1)
w4, = plt.plot(t, cross4)
label = fm.format(' {0} - {1} - Crosscorrelation',station2, station3)
plt.title(label)

plt.ylabel("Correlation coefficient")
plt.title(label)
plt.subplot(428, sharex=ax1)
maxvalue[0] = wf4.max()
sumvalue = sum(wf4)/samprate
for j in range (len(wf4)):
	if(maxvalue[0]-epsilon <= wf4[j]):
		maxtime[0] = float(j);
for j in range (len(wf4)):
	wf4[j] = wf4[j]/sumvalue
maxtime[0] = maxtime[0]/samprate
maxvalue[0] = maxvalue[0]/sumvalue


lapx0, lapstdv , sumd = gy.FitLaplace(wf4,samprate, halfwidth)
maxtime[1] = lapx0


for j in range(len(twf)):
	arg = 0.5*(1./lapstdv)*math.exp(-1.*abs(j/samprate-lapx0)/lapstdv)
	lapdist[j] = sumd*arg
	arg = (1./math.sqrt(2.*math.pi))*(1./lapstdv)*math.exp(-1.*(j/samprate-lapx0)*(j/samprate-lapx0)/2./lapstdv/lapstdv)
	gaudist[j] = sumd*arg


maxvalue[1] = max(lapdist)
plt.plot(twf, wf4, 'blue', maxtime, maxvalue,'ro', twf, lapdist, 'red', twf, gaudist, 'green')

label = fm.format(' {0} - {1} - CC Envelope',station2, station3)
plt.title(label)

label = fm.format('Seconds after {0}',start_time)
plt.xlabel(label)
plt.ylabel("Correlation coefficient")
plt.subplots_adjust(bottom=0.2)

figurename = triad+'_'+start_time+'.png'
pylab.savefig(figurename)
plt.show()
