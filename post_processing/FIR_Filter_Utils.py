#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 14:52:52 2018

@author: jmotyka
"""

from scipy import signal
import numpy as np
import math as math

#from scipy.fftpack import fft, fftshift
#import matplotlib.pyplot as plt

def SymmetricFIR(N,winType,rawData):
    if winType == 'triang':
        tmp = signal.triang(N,True)
        w = np.divide(tmp, np.sum(tmp))

    numPts = len(rawData)

    filtData = rawData
    for k in range(numPts):
        if k < int(N/2):
            filtData[k] = rawData[k]
            #print('low' )
        elif k >= (numPts - math.floor(N/2)):
            filtData[k] = rawData[k]
            #print('high' )
        else:
            tmp = rawData[k-math.floor(N/2):k+math.ceil(N/2)]
            filtData[k] = np.dot(w,tmp)

    return filtData


