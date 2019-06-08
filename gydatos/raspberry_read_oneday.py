# -*- coding: utf-8 -*-
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import math


starttime=UTCDateTime.utcnow()
#
# Get epoch time from UTCDateTime
# 
day_length = 24.*3600.
epochtime = starttime.timestamp
print(epochtime)
#
# Set start time at beginning of previous day
#
numdays=epochtime/day_length-1
print(numdays)
start_epoch=day_length*math.floor(numdays)
print(start_epoch)
starttime=UTCDateTime(start_epoch)
strpnow=UTCDateTime.strftime(starttime,'%Y %j')
print(strpnow)
#
# Get one day of data
#
endtime=starttime+day_length
client = Client(base_url='https://fdsnws.raspberryshakedata.com/')
waveform0 = client.get_waveforms('AM','R8A6C','00','EHZ',starttime,endtime)

n = str(waveform0)
print('Number of traces:',n,"\n")
for tr in waveform0:
     tr.data = np.require(tr.data, dtype=np.int32)
waveform0.write('/home/ShareData/'+strpnow.replace(" ","")+'.mseed',format='MSEED')
fname='OneDay'
startplot=starttime
waveform0.plot(starttime=startplot,outfile=fname+'_0.png')



