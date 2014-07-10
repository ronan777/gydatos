#!/usr/bin/env Python

# This module contains a number of general-use utilities 
# to read CSS data and do limited, and extendable, signal 
# processing.
#    
# 
# Developed June 2014 - Ronan Le Bras - Gydatos LLC
#

# import all necessary modules

import struct
import math
import array as ar
import matplotlib.pyplot as plt
from   matplotlib.widgets import Button, Slider
import string
import time
import calendar
#import scipy as sc
#import scipy.sparse
#import numpy as np
#import scipy.io.matlab.byteordercodes
from pylab import *
#import pyaudio
import wave
import glob
import os
import collections
import subprocess
import matplotlib.mpl as mp
import matplotlib.cm as cm
import ConfigParser
import cx_Oracle
import pylab

fm = string.Formatter()

Site = collections.namedtuple('Site', ['sta','lat','lon'], verbose=True)
Detection = collections.namedtuple('Detection',
['sta','Det_num','time','dt_cc','duration','freq','whale_disc','pow_15_25','pow_30_50','pow_50_60','spectim'], 
verbose=True)
Association = collections.namedtuple('Association',
['sta','Det_num','Ev_num','time','delta','StaEvAz'],verbose=True)
InfraSta =('I02AR','I04AU','I05AU','I06AU','I07AU','I08BO','I09BR','I10CA','I11CV','I13CL','I14CL','I17CI','I18DK','I19DJ','I21FR','I22FR','I23FR','I24FR','I26DE','I27DE','I30JP','I31KZ','I32KE','I33MG','I34MN','I35NA','I36NZ','I37NO','I39PW','I41PY','I42PT','I43RU','I44RU','I45RU','I46RU','I47ZA','I48TN','I49GB','I50GB','I51GB','I52GB','I53US','I55US','I56US','I57US','I58US','I59US')
#InfraSta = ('I52GB','I53US','I55US','I56US','I57US','I58US','I59US')
#InfraSta = ('I02AR','I04AU')

KMperDEG = 111.32
LARGE = 9999999.
epsilon = 0.001

config = ConfigParser.RawConfigParser()
config.read('/home/smd/lebras/polarinfra.cfg')

station   = config.get('Polar','station')
account = config.get('Polar', 'account')


con = cx_Oracle.connect(account)
print ('Version:',con.version)

cur = con.cursor()

for i in range(len(InfraSta)):
	fig = plt.figure(figsize=(12,12))
	querysel3valid = fm.format(' Select sta,arid,orid,delta,seaz from sel3.assoc@adb where delta between 20. and 60. and sta=\'{0}\' and arid not in (select arid from validinfra)', InfraSta[i])
	cur.execute(querysel3valid)
	N = 0
	radii = ar.array('f',[])
	az = ar.array('f',[])
	
	for result in cur:
		N = N+1
		print (N, result)
		sta,arid,orid,delta,seaz = result
		radii.append(float(delta))
		if(seaz > 0. and seaz < 90.):
			azz = 90.-seaz
		if(seaz > 90. and seaz < 360.):
			azz = 450.-seaz
		az.append(2.*pi*float(azz)/360.)
	if(N == 0):
		continue	
	area=20.
	colors=az	
	ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='polar')
	c = plt.scatter(az, radii, color='grey', s=area, cmap=cm.hsv)
	querysel3nonvalid = fm.format(' Select sta,arid,orid,delta,seaz from sel3.assoc@adb where delta between 20. and 60. and sta=\'{0}\' and arid in (select arid from validinfra)', InfraSta[i])
	cur.execute(querysel3nonvalid)
	N = 0
	radii = ar.array('f',[])
	az = ar.array('f',[])
	
	for result in cur:
		N = N+1
		print (N, result)
		sta,arid,orid,delta,seaz = result
		radii.append(float(delta))
		if(seaz > 0. and seaz < 90.):
			azz = 90.-seaz
		if(seaz > 90. and seaz < 360.):
			azz = 450.-seaz
		az.append(2.*pi*float(azz)/360.)
	if(N == 0):
		continue	
	area=20.
	colors=az	
	ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='polar')
	c = plt.scatter(az, radii, color='red', s=area, cmap=cm.hsv)
	c.set_alpha(0.75)

	#ax.set_theta_direction(-1)
	label = fm.format('SEL3 Station {0} ',sta)	
	plt.xlabel(label)
	
	figurename = fm.format('SEL3Polar{0}',sta)
	pylab.savefig(figurename, bbox_inches='tight')
#	show()

cur.close()
con.close()

