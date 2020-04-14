#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 15:52:17 2018

@author: jmotyka
"""

import numpy as np
import math as math


# Create some constants that specify array indices
DMU381  = int(0)
NOVATEL = int(1)
RTK     = int(2)

# Define constants used in the function
X_AXIS      = int(0)
Y_AXIS      = int(1)
Z_AXIS      = int(2)
NUM_OF_AXES = int(3)

Q0 = int(0)
Q1 = int(1)
Q2 = int(2)
Q3 = int(3)

ROLL  = int(0)
PITCH = int(1)
YAW   = int(2)

LAT = int(0)
LON = int(1)
ALT = int(2)

N_AXIS = int(0)
E_AXIS = int(1)
D_AXIS = int(2)

# Calibration constants
SLOPE  = int(0)
OFFSET = int(1)
AMP    = int(2)
RISE   = int(3)
SHIFT  = int(4)

# Set the median and standard deviation (based on nominal accelerometer
#   calibration values)
MU    = 4190000
SIGMA = 300000

# General conversion constants
DegToRad = ( np.pi/180.0 )
RadToDeg = ( 180.0/np.pi )

GRAV_ACCEL = 9.80655

FEET_TO_METERS = 0.3048

MILLI_G_TO_G = 1.0e-3

TWO_PI = 2.0*math.pi

