#!/usr/bin/env Python

# This module test fft for calibration of amplitude


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


fft1 = abs(scipy.fft(wf1))
freqs = scipy.fftpack.fftfreq(len(wf1),1./samprate)		
label = 'Spectrum '+pr_time
print(label)
plt.title(label)

plt.plot(time1,power2,'ro-')
	
plt.show()

