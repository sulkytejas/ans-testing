import os
import sys

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import numpy as np
import copy
import math

from PeakDetection import *
from Constants import *

D2R = math.pi/180.0

class sysEvaluation:
    ''' compare the performances of two INS systems: test, reference
    '''
    def __init__(self):
        ''' initialize the class
        '''
        self.leapseconds        = 18  #2019
        self.year               = 0
        self.month              = 0
        self.day                = 0
        self.sow                = 0
        
        # menu input
        self.ReferenceFileName  = ''
        self.TestFileName       = ''
        self.timeStr            = ''
        self.dirname            = ''
        self.bErrorInReading    = False

        # output filenames
        self.trajectoryFileName = ''
        self.attitudesFileName  = ''
        self.attDifFileName     = ''
        self.filesuffix         = ''        

        self.RefKMLFileName     = ''
        self.TestKMLFileName    = ''
        self.bGeneratePDF       = True
        self.RefKMLFileName     = ''
        self.TestKMLFileName    = ''
        self.bTestGPSavailable  = False
        self.bTestHeadingAvailable = False

        self.ReferenceSensor    = ''
        self.TestSensor         = ''
        self.bCPT7asReference   = False
        self.bUseBestGNSSPos    = False
        
        # week, gpssecond, latitude, longitude, altitude, roll, pitch, azimuth, N, E, H
        self.ReferenceData      = []
        self.TestData           = []
        self.ReferenceRate      = 20  # Hz
        self.magtestData        = []

        # ECEF coordinates and rotation matrix R_E^N of the start point
        self.startECEF          = []
        self.startRNE           = []
        self.firstRec           = True

        # time and initial attitude alignment, need to change for different data set
        self.bAutoSynAndAlign   = False
        
        self.ReferenceTimeAlign =  0
        self.TestTimeAlign      =  0
        self.initialTimeForRollPitch   = 0
        self.initialTimeForYaw         = 0
        self.interestBegin      = 0.0
        self.interestEnd        = 999999.0

        self.Differences        = []  
        self.strRefFixedRate    = ''
        self.strTestFixedRate   = ''


    def showmedifferences(self):
        if not self.Differences:
            if self.bAutoSynAndAlign == True:
                self.autoSynchronize()
            else:
                self.timeAlignment()
            self.initialAttitudesAlignment()
            self.solutionDifferences()
        self.TrajectoryPlot()
        return

    def showmeAttitudes(self):
        if not self.Differences.any():
            if self.bAutoSynAndAlign == True:
                self.autoSynchronize()
            else:
                self.timeAlignment()
            self.initialAttitudesAlignment()
            self.solutionDifferences()
        self.AttitudesPlot()
        return

    def showmeAttitudesDifferences(self):
        if not self.Differences.any():
            if self.bAutoSynAndAlign == True:
                self.autoSynchronize()
            else:
                self.timeAlignment()
            self.initialAttitudesAlignment()
            self.solutionDifferences()
        self.AttitudesDifPlot()
        return

    def timeAlignment(self):
        ''' openIMU solution does not have GPS time tag
        '''
        self.TestData[:,1] += self.ReferenceTimeAlign - self.TestTimeAlign

    def autoSynchronize(self):
        ''' synchronize two series automatically
        '''
        sig = [SignalPeaks() for _ in range(2)]

        sig[0].set_name('Novatel')
        sig[0].set_signal(self.ReferenceData[:,1], np.abs(self.ReferenceData[:,6]), 20)
        sig[0].set_smoothingParameters(40,'triang')
        sig[0].set_peakSeparationTime(0.25)
        sig[0].set_widthThreshold(12.5)
        sig[0].set_amplitudeThreshold(6.0)
        sig[0].SmoothData()
        sig[0].FindPeaks()
        sig[0].FindSeparationTimes()

        sig[1].set_name('Test')
        sig[1].set_signal(self.TestData[:,1], np.abs(self.TestData[:,6]),100)
        sig[1].set_smoothingParameters(40,'triang')
        sig[1].set_peakSeparationTime(0.25)
        sig[1].set_widthThreshold(12.5)
        sig[1].set_amplitudeThreshold(6.0)
        sig[1].SmoothData()
        sig[1].FindPeaks()
        sig[1].FindSeparationTimes()

        deltaT = sig[0].FindOptimalAlignment(sig[1])
        dT = deltaT[0,0]
        print('Time difference: '+str(dT))

        self.TestData[:,1] += dT

        
    def initialAttitudesAlignment(self):
        ''' align the intial attitudes for comparison between systems
        '''
        refRPIndex   = -1
        testRPIndex  = -1
        refYawIndex  = -1
        testYawIndex = -1
        if self.initialTimeForRollPitch   > 0.0:
            referencetime = np.fabs(self.ReferenceData[:,1] - self.initialTimeForRollPitch)
            index         = np.argmin(referencetime)
            if referencetime[index]<0.01:
                refRPIndex = index
            print('reference ', index,self.ReferenceData[index,1])
        
            tesTime = np.fabs(self.TestData[:,1] - self.initialTimeForRollPitch)
            index   = np.argmin(tesTime)
            if tesTime[index]<0.01:
                testRPIndex = index
            print('test data ', index,self.TestData[index,1])

        if self.bTestHeadingAvailable:
            if self.initialTimeForYaw != self.initialTimeForRollPitch:

                referencetime = np.fabs(self.ReferenceData[:,1] - self.initialTimeForYaw)
                index         = np.argmin(referencetime)
                if referencetime[index]<0.01:
                    refYawIndex = index
                print('reference ', index,self.ReferenceData[index,1])
            
                tesTime = np.fabs(self.TestData[:,1] - self.initialTimeForYaw)
                index   = np.argmin(tesTime)
                if tesTime[index]<0.01:
                    testYawIndex = index
                print('test data ', index,self.TestData[index,1])
            else:
                refYawIndex  = refRPIndex
                testYawIndex = testRPIndex

        if (refRPIndex==-1 or testRPIndex==-1 or refYawIndex==-1 or testYawIndex==-1):
            pass
        else:
            # week, gpssecond, latitude, longitude, altitude, roll, pitch, azimuth, N, E, H
            initialRolldif  = self.subangle(self.ReferenceData[refRPIndex][5], self.TestData[testRPIndex][5])
            initialPitchdif = self.subangle(self.ReferenceData[refRPIndex][6], self.TestData[testRPIndex][6])
            self.TestData[:,5] += initialRolldif
            self.TestData[:,6] += initialPitchdif 
            if self.bTestHeadingAvailable:
                initialYawdif   = self.subangle(self.ReferenceData[refYawIndex][7], self.TestData[testYawIndex][7])
                self.TestData[:,7] += initialYawdif
                normalIndex = (self.TestData[:,7] < -180.0)
                self.TestData[normalIndex,7] += 360.0
                normalIndex = (self.TestData[:,7] >= 180.0)    
                self.TestData[normalIndex,7] -= 360.0
         
        return

    def solutionDifferences(self):
        ''' compute the differences between two solutions
        '''
        result = []

        if self.ReferenceRate == 20:
            fiveIndex   = (np.floor(self.TestData[:,1]*100+0.5)%5 == 0) # to match 20 hz reference data
            newtestData = self.TestData[fiveIndex,:]
        else:
            newtestData = self.TestData

        # week, gpssecond, latitude, longitude, altitude, roll, pitch, azimuth, N, E, H
        iTestIndex    = 0
        timeThreshold = 0.001
        for j in range(len(self.ReferenceData)):

            for i in range(iTestIndex+1,len(newtestData)):
                dt = newtestData[i][1] - self.ReferenceData[j][1]
                if math.fabs(dt) < timeThreshold:
                    epoch    = self.ReferenceData[j][1]
                    rolldif  = self.subangle(newtestData[i][5], self.ReferenceData[j][5])
                    pitchdif = self.subangle(newtestData[i][6], self.ReferenceData[j][6])
                    yawdif   = 0.0
                    if self.bTestHeadingAvailable:
                        yawdif = self.subangle(newtestData[i][7], self.ReferenceData[j][7])
                    nDif     = newtestData[i][8] - self.ReferenceData[j][8]
                    eDif     = newtestData[i][9] - self.ReferenceData[j][9]
                    dDif     = newtestData[i][10] - self.ReferenceData[j][10]
                    hDis     = nDif*nDif+eDif*eDif
                    hvDis    = math.sqrt(hDis+dDif*dDif)
                    hDis     = math.sqrt(hDis)
                    temp     = [epoch, rolldif, pitchdif, yawdif, nDif, eDif, dDif, hDis, hvDis]
                    result.append(temp)
                    iTestIndex = i
                    break
                if (iTestIndex != 0 and dt > 2.0):
                    break
                
        self.Differences = []
        self.Differences = copy.deepcopy(np.array(result))

        return
        
    def readAceinnaFile(self):
        ''' read Aceinna openIMU or openIMU+GPS solution
        '''
        insFileName = self.TestFileName
        bIdPacket   = False
        bA1Packet   = False
        bE2Packet   = False
        b300RI_VGAH = False
        bOpenIMU335 = False
        self.strTestFixedRate   = ''

        try:
            with open(insFileName) as afile:
                line = afile.readline()
                #print(line)
                if line.find('gps_update, gps_valid')>=0:
                    self.bTestGPSavailable     = True
                    self.TestSensor            = 'Ublox_300ZI_e2'
                    self.bTestHeadingAvailable = True
                    bIdPacket                  = True
                elif line.find('timeCntr,time,roll,pitch,heading,xAccel,yAccel,zAccel')>=0:
                    self.bTestGPSavailable     = True
                    self.TestSensor            = 'Ublox_300ZI_e2'
                    self.bTestHeadingAvailable = True
                    bE2Packet                  = True
                elif line.find('time,roll,pitch,xRate,yRate,zRate,xAccel,yAccel,zAccel,!!!')>=0:
                    self.bTestGPSavailable     = False
                    self.TestSensor            = '300RI-VG_AHRS'
                    self.bTestHeadingAvailable = False
                    b300RI_VGAH                = True
                elif line.find('time,roll,pitch,xRate,yRate,zRate,xAccel,yAccel,zAccel,lat,long,gps_yaw,gps_roll,gps_pitch,cog,sog')>=0 or \
                     line.find('time,rsvd1,rsvd2,xRate,yRate,zRate,xAccel,yAccel,zAccel,lat,long,gps_yaw,gps_pitch,gps_roll,cog,sog')>=0:
                    self.bTestGPSavailable     = False
                    self.TestSensor            = '300RI_VG-AHRS_INS'
                    self.bTestHeadingAvailable = True
                    b300RI_VGAH                = False
                    b300RI_VGAH_INS            = True
                    print('b300RI_VGAH_INS')
                    
                elif line.find('counter, pitch, roll, yaw, rateX, rateY, rateZ, accelX, accelY, accelZ, magX, magY, magZ, temp, itow, BIT')>=0:
                    self.bTestGPSavailable     = False

                    self.bTestHeadingAvailable = True
                    bOpenIMU335                = True
                    if (insFileName.find('OpenIMU')>=0 or insFileName.find('openIMU')>=0 or insFileName.find('openimu')>=0):
                        self.TestSensor        = 'OpenIMU335'
                    else:
                        self.TestSensor        = 'MTLT'

                    if insFileName.find('mtlt')>=0:
                        self.TestSensor        = 'MTLT'
                    else:
                         self.TestSensor       = 'OpenIMU335'
                         
                else:
                    self.bTestGPSavailable = False
                    if line.find('heading')>=0:
                        self.TestSensor            = '300ZI_VG-AHRS'
                        self.bTestHeadingAvailable = True
                    elif line.find('roll,pitch,yaw,xRate,yRate,zRate,xAccel,yAccel,zAccel,xMag,yMag,zMag,xRateTemp,gpsITOW,BIT')>=0:
                        bA1Packet = True
                        basename  = os.path.basename(insFileName)
                        self.TestSensor            = basename[7:19]
                        self.bTestHeadingAvailable = False
                    else:
                        self.TestSensor            = 'OpenIMU300ZI_INS'
                        self.bTestHeadingAvailable = False
                        
                lines  = afile.readlines()
                result = []
                i      = 0
                for x in lines:
                    if len(x)>400: continue
                    try:
                        values = x.split(',')
                    except ValueError:
                        continue
                    else:
                        ned    = [0,0,0]
                        if bIdPacket:
                            sow    = float(values[1])
                            lat    = float(values[8])
                            lon    = float(values[9])
                            alt    = float(values[10])
                            roll   = float(values[14])
                            pitch  = float(values[15])
                            yaw    = self.normalangle(float(values[16]))                            
                            if math.fabs(lat)>0 and math.fabs(lon)>0:
                                LLA    = [lat*D2R, lon*D2R, alt]
                                ned, _ = self.LLA2Base( LLA, self.startECEF, self.startRNE )
                        elif bE2Packet:
                            sow    = float(values[1])
                            roll   = float(values[2])
                            pitch  = float(values[3])
                            yaw    = float(values[4])
                            lat    = float(values[23])
                            lon    = float(values[24])
                            alt    = float(values[25])     
                        elif bA1Packet:
                            sow    = i/100.0
                            roll   = float(values[0])
                            pitch  = float(values[1])                            
                            yaw    = 0.0
                            lat    = 0
                            lon    = 0
                            alt    = 0
                            i      = i+1
                        elif bOpenIMU335:
                            sow    = float(values[0])/1000.0
                            roll   = float(values[1])
                            pitch  = float(values[2])                            
                            yaw    = float(values[3])                            
                            lat    = 0
                            lon    = 0
                            alt    = 0
                        elif b300RI_VGAH:
                            if math.fabs(float(values[1])) > 0.0 or math.fabs(float(values[2])) > 0.0:
                                sow   = float(values[0])-1575832700.0 # UNIX time
                                roll  = -float(values[1])
                                pitch = -float(values[2])
                                yaw   = 0.0
                                lat   = 0.0
                                lon   = 0.0
                                alt   = 0.0
                            else:
                                continue
                        elif b300RI_VGAH_INS:
                            if math.fabs(float(values[11])) > 0.0 or math.fabs(float(values[12])) > 0.0 or math.fabs(float(values[13])) > 0.0:
                                sow   = float(values[0])-1576300700.0
                                roll  = float(values[13])/D2R
                                pitch = float(values[12])/D2R
                                yaw   = self.normalangle(float(values[11])/D2R)
                                lat   = 0.0
                                lon   = 0.0
                                alt   = 0.0
                            else:
                                continue
                        else:
                            sow    = float(values[1])
                            roll   = float(values[2])
                            pitch  = float(values[3])
                            yaw    = 0.0
                            if self.bTestHeadingAvailable:
                                yaw = self.normalangle(float(values[4]))
                            lat    = 0
                            lon    = 0
                            alt    = 0
                            
                        pass
                    

                        if self.bTestGPSavailable == True:
                            if math.fabs(lat)>0 and math.fabs(lon)>0:
                                LLA    = [lat*D2R, lon*D2R, alt]
                                ned, _ = self.LLA2Base( LLA, self.startECEF, self.startRNE )

                        temp   = [0.0, sow, lat, lon, alt, roll, pitch, yaw, ned[0], ned[1], ned[2]]
                        if self.bTestGPSavailable == True:
                            if (sow>=self.interestBegin and sow<=self.interestEnd):
                                result.append(temp)
                        else:
                            result.append(temp)
                        
                self.TestData = []
                self.TestData = copy.deepcopy(np.array(result))

        except:
            pass

        return

    def readInsFile(self,sensor):
        ''' read in the information in INS solution file
            week, sow, lat, lon, roll, pitch, azimuth
        '''   
        insFileName = ''
        if sensor  == 'Novatel_FlexPak6_iMAR':
            insFileName = self.ReferenceFileName        
        elif sensor == 'Frii':
            insFileName = self.TestFileName
        elif sensor=='Novatel_CPT7':
            if self.bCPT7asReference:
                insFileName = self.ReferenceFileName
            else:
                insFileName = self.TestFileName

        bAfterFirstGood = False
        result          = []
        lines           = []
        
        if self.bErrorInReading==True:
            with open(insFileName) as afile:
                for cnt, line in enumerate(afile):
                    print("Line {}: {}".format(cnt+1,line))

        try:                        
            with open(insFileName) as afile:
                lines    = afile.readlines()
                nrecords = len(lines)
                
                print('Total lines: ',nrecords)
                
                i        = 0
                nGNSSpos = 0
                nGNSSfix = 0
                while i < nrecords:
                    x = lines[i]
                    i = i+1
                        
                    if bAfterFirstGood==True or x.find('INS_SOLUTION_GOOD')>=0:
                        bAfterFirstGood = True

                        try:
                            values = x.split(',')
                        except ValueError:
                            continue
                        else:
                            if (self.bUseBestGNSSPos == False and x.find('INSPVAXA')>=0):
                                week   = float(values[5] )
                                sow    = float(values[6] )
                                lat    = float(values[11])
                                lon    = float(values[12])
                                alt    = float(values[13]) + float(values[14])
                                roll   = float(values[18])
                                pitch  = float(values[19])
                                yaw    = float(values[20])

                            elif x.find('BESTGNSSPOSA')>=0:
                                nGNSSpos += 1
                                if self.bUseBestGNSSPos==True:
                                    week   = float(values[5] )
                                    sow    = float(values[6] )
                                    lat    = float(values[11])
                                    lon    = float(values[12])
                                    alt    = float(values[13]) + float(values[14])
                                    roll   = 0.0
                                    pitch  = 0.0
                                    yaw    = 0.0
                                if (x.find('NARROW_INT')>=0 and float(values[6])>=self.interestBegin and float(values[6])<=self.interestEnd):
                                    nGNSSfix += 1
                                if self.bUseBestGNSSPos==False: continue
                            elif (x.find('INSPVASA')>=0):
                                week   = float(values[1] )
                                sow    = float(values[3] )
                                lat    = float(values[4])
                                lon    = float(values[5])
                                alt    = float(values[6])
                                roll   = float(values[10])
                                pitch  = float(values[11])
                                yaw    = float(values[12])

                            else:
                                continue
                            
                            yaw        = self.normalangle(yaw)

                            if self.year == 0:
                                calendartime = self.timefromGPS(week,sow)
                                self.year    = calendartime[0]
                                self.month   = calendartime[1]
                                self.day     = calendartime[2]
                                self.sow     = sow
                                
                            if (sow>=self.interestBegin and sow<=self.interestEnd):
                                LLA = [lat*D2R, lon*D2R, alt]
                                if self.firstRec:
                                    self.startRNE  = np.array(self.RneFromLLA( LLA ))
                                    self.startECEF = np.array(self.LLA2ECEF( LLA ))
                                    ned = np.zeros( [NUM_OF_AXES,1] )
                                    self.firstRec = False
                                else:
                                    ned, _ = self.LLA2Base( LLA, self.startECEF, self.startRNE )
                                    
                                temp   = [week, sow, lat, lon, alt, roll, pitch, yaw, float(ned[0]), float(ned[1]), float(ned[2])]
                                #print(temp)
                                result.append(temp)

                if nGNSSpos > 1:
                    fixedRate  = (nGNSSfix/nGNSSpos)*100.0
                else:
                    fixedRate = 0
                if sensor=='Novatel_FlexPak6_iMAR':
                    self.ReferenceData = []
                    self.ReferenceData = copy.deepcopy(np.array(result))
                    interval = self.ReferenceData[1][1] - self.ReferenceData[0][1]
                    self.ReferenceRate = math.floor(1.0/interval + 0.5)
                    if (fixedRate>0.0):
                        self.strRefFixedRate  = 'RTK fix rate: ' + '{:.{prec}f}'.format(fixedRate, prec=1) + '%'
                        print(self.strRefFixedRate)
                elif (sensor=='Frii'):
                    self.TestData = []
                    self.TestData = copy.deepcopy(np.array(result))
                    if (fixedRate>0.0):
                        self.strTestFixedRate  = 'RTK fix rate: ' + '{:.{prec}f}'.format(fixedRate, prec=1) + '%'
                        print(self.strTestFixedRate)
                elif sensor == 'Novatel_CPT7':
                    if self.bCPT7asReference:
                        self.ReferenceData = []
                        self.ReferenceData = copy.deepcopy(np.array(result))
                        if (fixedRate>0.0):
                            self.strRefFixedRate  = 'RTK fix rate: ' + '{:.{prec}f}'.format(fixedRate, prec=1) + '%'
                            print(self.strRefFixedRate)
                    else:
                        self.TestData = []
                        self.TestData = copy.deepcopy(np.array(result))
                        if (fixedRate>0.0):
                            self.strTestFixedRate  = 'RTK fix rate: ' + '{:.{prec}f}'.format(fixedRate, prec=1) + '%'
                            print(self.strTestFixedRate)
        except:
            pass
        return
         
    def generateKML(self):
        self.writeKML(True)
        if  self.bTestGPSavailable:
            self.writeKML(False)
        

    def writeKML(self, bReference=True):      
        '''Writing the kml file
        '''
        fname = ''
        data  = []
        if bReference:
            fname = self.RefKMLFileName
            data  = self.ReferenceData
            color = 'ffff00'
        else:
            fname = self.TestKMLFileName
            data  = self.TestData
            color = '0000ff'

        with open(fname,'w') as outf:
            outf.write("<?xml version='1.0' encoding='UTF-8'?>\n")
            outf.write("<kml xmlns='http://www.opengis.net/kml/2.2'>\n")
            outf.write("<color>"+color+"</color>\n")
            outf.write("<scale>0.100</scale>\n")
            outf.write("<Icon>\n")
            outf.write("<href>http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png</href>\n")
            outf.write("</Icon>\n")
            outf.write("<Document>\n")
            i = 0
            for row in data:
                if (True or i%10 == 0):
                    outf.write("   <Placemark>\n")
                    outf.write("       <Point>\n")
                    outf.write("           <coordinates>" + "{:.{prec}f}".format(row[3],prec=11) + "," + "{:.{prec}f}".format(row[2],prec=11) + "," + "{:.{prec}f}".format(row[4],prec=4) + "</coordinates>\n")
                    outf.write("       </Point>\n")
                    outf.write("   </Placemark>\n")
                i += 1
            outf.write("</Document>\n")
            outf.write("</kml>\n")
            outf.close()
            

    def LLA2ECEF( self, LLA ):
        ''' convert lat, lon, and alititude to ECEF coordinates
        '''
        E_MAJOR  = 6.378137e+6
        f        = 1.0/298.257223563
        E_ECC_SQ = 2*f - f*f
        E_MINOR_OVER_MAJOR_SQ = 1 - E_ECC_SQ

        cosLat = math.cos( LLA[LAT] )
        sinLat = math.sin( LLA[LAT] )
        cosLon = math.cos( LLA[LON] )
        sinLon = math.sin( LLA[LON] )

        N = E_MAJOR / math.sqrt( 1.0 - ( E_ECC_SQ * sinLat * sinLat ) )

        ECEF = np.zeros([NUM_OF_AXES,1])
        ECEF[X_AXIS] = ( N + LLA[ALT] ) * cosLat * cosLon
        ECEF[Y_AXIS] = ( N + LLA[ALT] ) * cosLat * sinLon
        ECEF[Z_AXIS] = ( ( E_MINOR_OVER_MAJOR_SQ * N ) + LLA[ALT] ) * sinLat
        
        return ECEF

    def ECEF2LLA(self, ECEF):
        ''' convert ECEF coordinates to lat, lon in deegrees and altitude in meters
        '''
        a   	= 6378137.0
        f       = 1.0/298.257223563
        e2      = 2*f-f*f
        x       = ECEF[0]
        y       = ECEF[1]
        z       = ECEF[2]
        R       = math.sqrt(x*x + y*y + z*z)
        ang     = math.atan(math.fabs(z/math.sqrt(x*x+y*y)))
        if z<0.0:
            ang = -ang
        lat1    = ang
        Rw      = math.sqrt(1-e2*math.sin(lat1) * math.sin(lat1))
        lat     = math.atan(math.fabs(math.tan(ang)*(1+(a*e2*math.sin(lat1))/(z*Rw))))
        if  z<0.0:
            lat = -lat
        while math.fabs(lat-lat1)>1e-12:
            lat1 = lat
            Rw   = math.sqrt(1-e2*math.sin(lat1)*math.sin(lat1))
            lat  = math.atan(math.fabs(math.tan(ang)*(1+(a*e2*math.sin(lat1))/(z*Rw))))
            if (z<0.0):
                lat =-lat
        if lat>math.pi:
            lat = lat-2.0*math.pi
        if math.fabs(x)<1e-12:
            if y>=0.0:
                lon = math.pi/2.0
            else:
                lon = 3.0*math.pi/2.0
        else:
            lon = math.atan(math.fabs(y/x))
            if x>0.0:
                if y>=0.0:
                    lon = lon
                else:
                    lon = 2.0*math.pi-lon
            else:
                if y>=0.0:
                    lon = math.pi-lon
                else:
                    lon = math.pi+lon
        Rw      = math.sqrt(1-e2*math.sin(lat)*math.sin(lat))
        Rn      = a/Rw
        ht      = R*math.cos(ang)/math.cos(lat)-Rn
        if lon>math.pi:
            lon = lon-2.0*math.pi
        LLA     = np.zeros([NUM_OF_AXES,1])
        LLA[0]  = lat / D2R
        LLA[1]  = lon / D2R
        LLA[2]  = ht
        return LLA

        
    def RneFromLLA( self, LLA ):
        sinLat = math.sin( LLA[LAT] )
        cosLat = math.cos( LLA[LAT] )
        sinLon = math.sin( LLA[LON] )
        cosLon = math.cos( LLA[LON] )

        Rne = np.zeros( [NUM_OF_AXES,NUM_OF_AXES] )
        Rne[X_AXIS,X_AXIS] = -sinLat * cosLon
        Rne[X_AXIS,Y_AXIS] = -sinLat * sinLon
        Rne[X_AXIS,Z_AXIS] =  cosLat

        Rne[Y_AXIS,X_AXIS] = -sinLon
        Rne[Y_AXIS,Y_AXIS] =  cosLon
        Rne[Y_AXIS,Z_AXIS] =  0.0

        Rne[Z_AXIS,X_AXIS] = -cosLat * cosLon
        Rne[Z_AXIS,Y_AXIS] = -cosLat * sinLon
        Rne[Z_AXIS,Z_AXIS] = -sinLat

        return Rne

    def LLA2Base( self, LLA, BaseECEF, Rne ):
        ''' Compute the ECEF coordinates (in the ECEF-frame) given LLA
        '''
        newECEF = self.LLA2ECEF( LLA )

        delta = np.zeros([NUM_OF_AXES,1])
        delta[X_AXIS] = newECEF[X_AXIS] - BaseECEF[X_AXIS]
        delta[Y_AXIS] = newECEF[Y_AXIS] - BaseECEF[Y_AXIS]
        delta[Z_AXIS] = newECEF[Z_AXIS] - BaseECEF[Z_AXIS]
        
        NED = np.zeros( [NUM_OF_AXES,1] )
        NED[X_AXIS] = Rne[X_AXIS,X_AXIS] * delta[X_AXIS] + \
                      Rne[X_AXIS,Y_AXIS] * delta[Y_AXIS] + \
                      Rne[X_AXIS,Z_AXIS] * delta[Z_AXIS]
        NED[Y_AXIS] = Rne[Y_AXIS,X_AXIS] * delta[X_AXIS] + \
                      Rne[Y_AXIS,Y_AXIS] * delta[Y_AXIS] + \
                      Rne[Y_AXIS,Z_AXIS] * delta[Z_AXIS]
        NED[Z_AXIS] = Rne[Z_AXIS,X_AXIS] * delta[X_AXIS] + \
                      Rne[Z_AXIS,Y_AXIS] * delta[Y_AXIS] + \
                      Rne[Z_AXIS,Z_AXIS] * delta[Z_AXIS]
                      
        return NED, newECEF

    def Open_Ref_Test_Files(self, namestr):
        ''' open the menu, get filename, read file, and plot
        '''
        sensor   = ''
        basename = os.path.basename(namestr)
        if self.dirname == '':
            self.dirname = os.path.dirname(namestr)
            if len(self.dirname) != 0 and self.dirname[-1] != '/':
                self.dirname += '/'
            if self.timeStr == '':
                dasloc = basename.find('-')
                dotloc = basename.find('.')
                self.timeStr = basename[dasloc+1:dotloc]
        
        if basename.find('Novatel_FLX6') >= 0 or basename.find('novatel_FLX6') >= 0:
            sensor                 = 'Novatel_FlexPak6_iMAR'
            self.ReferenceFileName = namestr
            self.ReferenceSensor   = 'Novatel_FlexPak6_iMAR'
            if self.filesuffix == '':
                self.RefKMLFileName    = self.dirname+sensor+'-'+self.timeStr+'.kml'
            else:
                self.RefKMLFileName    = self.dirname+sensor+'-'+self.timeStr+'-'+self.filesuffix+'.kml'
        elif basename.find('novatel_CPT7') >= 0:

            sensor                         = 'Novatel_CPT7'
            if self.bCPT7asReference:
                self.ReferenceFileName     = namestr
                self.ReferenceSensor       = 'Novatel_CPT7'
                if self.filesuffix == '':
                    self.RefKMLFileName    = self.dirname+sensor+'-'+self.timeStr+'.kml'
                else:
                    self.RefKMLFileName    = self.dirname+sensor+'-'+self.timeStr+'-'+self.filesuffix+'.kml'
            else:
                self.TestFileName          = namestr
                self.TestSensor            = 'Novatel_CPT7'
                self.bTestGPSavailable     = True
                self.bTestHeadingAvailable = True                
                if self.filesuffix == '':
                    self.TestKMLFileName   = self.dirname+sensor+'-'+self.timeStr+'.kml'
                else:
                    self.TestKMLFileName   = self.dirname+sensor+'-'+self.timeStr+'-'+self.filesuffix+'.kml'
                    
        elif basename.find('Frii') >= 0 or basename.find('frii') >= 0 or basename.find('ins2000') >= 0:
            sensor                 = 'Frii'
            self.TestFileName      = namestr
            self.TestSensor        = 'INS2000'
            self.bTestGPSavailable     = True
            self.bTestHeadingAvailable = True 
            if self.filesuffix == '':
                self.TestKMLFileName   = self.dirname+sensor+'-'+self.timeStr+'.kml'
            else:
                self.TestKMLFileName   = self.dirname+sensor+'-'+self.timeStr+'-'+self.filesuffix+'.kml'
                
        elif basename.find('data') >= 0 or basename.find('log') >= 0 or basename.find('txt') >= 0:
            sensor                 = 'Aceinna'
            self.TestFileName      = namestr
            if self.filesuffix == '':
                self.TestKMLFileName   = self.dirname+sensor+'-'+self.timeStr+'.kml'
            else:
                self.TestKMLFileName   = self.dirname+sensor+'-'+self.timeStr+'-'+self.filesuffix+'.kml'               
            self.readAceinnaFile()

        if self.filesuffix == '':
            self.trajectoryFileName    = self.dirname+'Trajectory-'+self.ReferenceSensor+'-'+self.TestSensor+'-'+self.timeStr+'.pdf'
            self.attitudesFileName     = self.dirname+'Attitudes-' +self.ReferenceSensor+'-'+self.TestSensor+'-'+self.timeStr+'.pdf'
            self.attDifFileName        = self.dirname+'AttitudesDiff-' +self.ReferenceSensor+'-'+self.TestSensor+'-'+self.timeStr+'.pdf'
        else:
            self.trajectoryFileName    = self.dirname+'Trajectory-'+self.ReferenceSensor+'-'+self.TestSensor+'-'+self.timeStr+'-'+self.filesuffix+'.pdf'
            self.attitudesFileName     = self.dirname+'Attitudes-' +self.ReferenceSensor+'-'+self.TestSensor+'-'+self.timeStr+'-'+self.filesuffix+'.pdf'
            self.attDifFileName        = self.dirname+'AttitudesDiff-' +self.ReferenceSensor+'-'+self.TestSensor+'-'+self.timeStr+'-'+self.filesuffix+'.pdf'
            
        if sensor != 'Aceinna':
            self.readInsFile(sensor)

        if sensor == 'Novatel_FlexPak6_iMAR':
            outfilename = self.dirname+self.ReferenceSensor+'-'+self.timeStr+'.pdf'
            self.singleSystemPlot(self.ReferenceData,outfilename,True,'g', self.ReferenceSensor, True)
        elif (sensor == 'Frii' or sensor == 'Aceinna' or sensor == 'UbloxF9_Swift' or sensor == 'FRII_D_PLUS'):
            outfilename = self.dirname+self.TestSensor+'-'+self.timeStr+'.pdf'
            self.singleSystemPlot(self.TestData,outfilename,True,'b',self.TestSensor,False)
        elif (sensor == 'Novatel_CPT7'):
            if self.bCPT7asReference:
                outfilename = self.dirname+self.ReferenceSensor+'-'+self.timeStr+'.pdf'
                self.singleSystemPlot(self.ReferenceData,outfilename,True,'g', self.ReferenceSensor)
            else:
                outfilename = self.dirname+self.TestSensor+'-'+self.timeStr+'.pdf'
                self.singleSystemPlot(self.TestData,outfilename,True,'b',self.TestSensor,False)

        return 

    def TrajectoryPlot(self):
        ''' plot the trajectory of two systems and their position differeneces
            and write the result to a file
        '''
        # epoch, rolldif, pitchdif, yawdif, nDif, eDif, dDif, hDis, hvDis
        if self.bTestGPSavailable:
            accutime = self.Differences[:,0]-self.Differences[0][0]
            north = self.Differences[:,4]
            east  = self.Differences[:,5]
            down  = self.Differences[:,6]
            length3D = self.Differences[:,8]
            meanN = np.mean(north)
            stdN  = np.std(north)
            meanE = np.mean(east)
            stdE  = np.std(east)
            meanD = np.mean(down)
            stdD  = np.std(down)
            meanL = np.mean(length3D)
            stdL  = np.std(length3D)

            strN  = 'mean={:.{prec}f}'.format(meanN,prec=2)+', std={:.{prec}f}'.format(stdN,prec=2)
            strE  = 'mean={:.{prec}f}'.format(meanE,prec=2)+', std={:.{prec}f}'.format(stdE,prec=2)
            strD  = 'mean={:.{prec}f}'.format(meanD,prec=2)+', std={:.{prec}f}'.format(stdD,prec=2)
            strL  = 'mean={:.{prec}f}'.format(meanL,prec=2)+', std={:.{prec}f}'.format(stdL,prec=2)

        # week, gpssecond, latitude, longitude, altitude, roll, pitch, azimuth, N, E, H
        with PdfPages(self.trajectoryFileName) as pdf:
            plt.figure(figsize=(7, 10))
            plt.plot(self.ReferenceData[:,9],self.ReferenceData[:,8], marker=',', linestyle='None', color='g', markersize=3,label=self.ReferenceSensor, rasterized=True)
            if self.bTestGPSavailable:
                plt.plot(self.TestData[:,9],self.TestData[:,8], marker=',', linestyle='None', color='b', markersize=3,label=self.TestSensor, rasterized=True)
            plt.axis('equal')
            plt.xlabel('East [m]')
            plt.ylabel('North [m]')
            plt.title('Trajectory')
            plt.legend()
            if self.bGeneratePDF == True:
                pdf.savefig()
            plt.close()

            plt.figure(figsize=(7, 10))
            ax = plt.subplot(3,1,1)
            plt.plot(self.ReferenceData[:,1],self.ReferenceData[:,8], marker=',', linestyle='None', color='g',markersize=3,label=self.ReferenceSensor, rasterized=True)
            plt.ylabel('[m]')
            plt.title('North')
            plt.legend()

            ax = plt.subplot(3,1,2)
            plt.plot(self.ReferenceData[:,1],self.ReferenceData[:,9], marker=',', linestyle='None', color='b',markersize=3,label=self.ReferenceSensor, rasterized=True)
            plt.ylabel('[m]')
            plt.title('East')
            plt.legend()

            ax = plt.subplot(3,1,3)
            plt.plot(self.ReferenceData[:,1],self.ReferenceData[:,10], marker=',', linestyle='None', color='r',markersize=3,label=self.ReferenceSensor, rasterized=True)
            plt.xlabel('Seconds of Week')
            plt.ylabel('[m]')
            plt.title('Up')
            plt.legend()
            ax.annotate(self.timeStr, xy=(-0.1, -0.4), xycoords='axes fraction',size=9)   
            if self.bGeneratePDF == True:
                pdf.savefig()
            plt.close()
                        
            if self.bTestGPSavailable:
                plt.figure(figsize=(7, 10))
                ax = plt.subplot(3,1,1)
                plt.plot(accutime,north, marker=',', linestyle='None', color='g', rasterized=True)
                plt.ylabel('[m]')
                plt.title('North Differences')
                ax.annotate(strN, xy=(0.70, 1.02), xycoords='axes fraction',size=10)            

                ax = plt.subplot(3,1,2)
                plt.plot(accutime,east, marker=',', linestyle='None', color='b', rasterized=True)
                plt.ylabel('[m]')
                plt.title('East Differences')
                ax.annotate(strE, xy=(0.70, 1.02), xycoords='axes fraction',size=10)

                ax = plt.subplot(3,1,3)
                plt.plot(accutime,down, marker=',', linestyle='None', color='r', rasterized=True)
                plt.xlabel('Seconds')
                plt.ylabel('[m]')
                plt.title('Up Differences')
                ax.annotate(strD, xy=(0.70, 1.02), xycoords='axes fraction',size=10)
                ax.annotate(self.timeStr, xy=(-0.1, -0.4), xycoords='axes fraction',size=9)   
                if self.bGeneratePDF == True:
                    pdf.savefig()
                plt.close()

                # 3D differences
                plt.figure(figsize=(7, 10))
                ax = plt.subplot(1,1,1)
                plt.plot(accutime,length3D, marker=',', linestyle='None', color='b', rasterized=True)
                plt.xlabel('Seconds')
                plt.ylabel('[m]')
                plt.title('3D Distances')
                ax.annotate(strL, xy=(0.70, 1.02), xycoords='axes fraction',size=10)
                ax.annotate(self.timeStr, xy=(-0.1, -0.4), xycoords='axes fraction',size=9)   
                if self.bGeneratePDF == True:
                    pdf.savefig()
                plt.close()  
                
                if self.bGeneratePDF == True:
                    pdf.savefig()
                plt.close()
            
        return

    def AttitudesPlot(self):
        ''' plot the attitudes of two systems
            and write the result to a file
        '''
        # week, gpssecond, latitude, longitude, altitude, roll, pitch, azimuth, N, E, H
        with PdfPages(self.attitudesFileName) as pdf:
            fig = plt.figure(figsize=(7, 10))
            plt.subplot(3,1,1)
            plt.plot(self.ReferenceData[:,1],self.ReferenceData[:,5], marker='.', linestyle='None', color='g', markersize=3,label=self.ReferenceSensor, rasterized=True)
            plt.plot(self.TestData[:,1],self.TestData[:,5], marker='.', linestyle='None', color='b', markersize=3,label=self.TestSensor, rasterized=True)
            plt.ylabel('[Deg]')
            plt.title('Roll')
            plt.legend()
            bx = plt.gca()
            bx.format_coord = lambda x,y: '%12f, %12f' % (x,y)

            ax = plt.subplot(3,1,2)
            plt.plot(self.ReferenceData[:,1],self.ReferenceData[:,6], marker='.', linestyle='None', color='g', markersize=3,label=self.ReferenceSensor, rasterized=True)
            plt.plot(self.TestData[:,1],self.TestData[:,6], marker='.', linestyle='None', color='b', markersize=3,label=self.TestSensor, rasterized=True)
            plt.ylabel('[Deg]')
            plt.title('Pitch')
            plt.legend()
            bx = plt.gca()
            bx.format_coord = lambda x,y: '%12f, %12f' % (x,y)
            if self.bTestHeadingAvailable == False:
                plt.xlabel('Seconds of Week')

            if self.bTestHeadingAvailable:
                ax = plt.subplot(3,1,3)
                plt.plot(self.ReferenceData[:,1],self.ReferenceData[:,7], marker='.', linestyle='None', color='g', markersize=3,label=self.ReferenceSensor, rasterized=True)
                plt.plot(self.TestData[:,1],self.TestData[:,7], marker='.', linestyle='None', color='b', markersize=3,label=self.TestSensor, rasterized=True)
                plt.xlabel('Seconds of Week')
                plt.ylabel('[Deg]')
                plt.title('Azimuth')
                plt.legend()
                bx = plt.gca()
                bx.format_coord = lambda x,y: '%12f, %12f' % (x,y)               
            ax.annotate(self.timeStr, xy=(-0.1, -0.4), xycoords='axes fraction',size=9) 
            plt.close()   
            if self.bGeneratePDF == True:
                pdf.savefig(fig)
        return

    

    def AttitudesDifPlot(self):
        ''' plot the attitude differences of two systems
            and write the result to a file
        '''
        # epoch, rolldif, pitchdif, yawdif, nDif, eDif, dDif
        accutime = self.Differences[:,0]-self.Differences[0][0]
        roll  = self.Differences[:,1]
        pitch = self.Differences[:,2]

        meanroll = np.mean(roll)
        stdroll  = np.std(roll)
        meanpitch = np.mean(pitch)
        stdpitch  = np.std(pitch)

        strroll   = 'mean={:.{prec}f}'.format(meanroll,prec=2)+', std={:.{prec}f}'.format(stdroll,prec=2)
        strpitch  = 'mean={:.{prec}f}'.format(meanpitch,prec=2)+', std={:.{prec}f}'.format(stdpitch,prec=2)

        if self.bTestHeadingAvailable:
            yaw   = self.Differences[:,3]
            meanyaw = np.mean(yaw)
            stdyaw  = np.std(yaw)
            stryaw    = 'mean={:.{prec}f}'.format(meanyaw,prec=2)+', std={:.{prec}f}'.format(stdyaw,prec=2)



        with PdfPages(self.attDifFileName) as pdf:
            fig = plt.figure(figsize=(7, 10))
            ax = plt.subplot(3,1,1)
            plt.plot(accutime,roll, marker=',', linestyle='None', color='g', markersize=3, rasterized=True)
            plt.ylabel('[Deg]')
            plt.title('Roll Differences')
            ax.annotate(strroll, xy=(0.70, 1.02), xycoords='axes fraction',size=10)

            ax = plt.subplot(3,1,2)
            plt.plot(accutime,pitch, marker=',', linestyle='None', color='b', markersize=3, rasterized=True)
            plt.ylabel('[Deg]')
            plt.title('Pitch Differences')
            ax.annotate(strpitch, xy=(0.70, 1.02), xycoords='axes fraction',size=10)
            if self.bTestHeadingAvailable == False:
                plt.xlabel('Seconds')


            if self.bTestHeadingAvailable:
                ax = plt.subplot(3,1,3)
                plt.plot(accutime,yaw, marker=',', linestyle='None', color='m', markersize=3, rasterized=True)
                plt.xlabel('Seconds')
                plt.ylabel('[Deg]')
                plt.title('Azimuth Differences')
                ax.annotate(stryaw, xy=(0.70, 1.02), xycoords='axes fraction',size=10)
                
            ax.annotate(self.timeStr, xy=(-0.1, -0.4), xycoords='axes fraction',size=9)
            plt.close()
            if self.bGeneratePDF == True:
                pdf.savefig(fig)
        return
    
    def singleSystemPlot(self,data,outputfilename,showGpsseconds,colorStyle,legendLabel, bReference=True):
        ''' plot the roll, pitch, and azimuth of both systems
        '''
        accutime = []
        if showGpsseconds == True:
            accutime = data[:,1]
        else:
            accutime = data[:,1]-data[0][1]
        wholeroll  = data[:,5]
        wholepitch = data[:,6]
        if bReference or self.bTestHeadingAvailable:
            wholeyaw   = data[:,7]

        with PdfPages(outputfilename) as pdf:
            fig = plt.figure(figsize=(7, 10))
            ax = plt.subplot(3,1,1)
            plt.plot(accutime,wholeroll, marker='.', linestyle='None', color=colorStyle, markersize=3,label=legendLabel, rasterized=True) 
            plt.ylabel('[Deg]')
            plt.title('Roll')
            plt.legend()
            if bReference:
                ax.annotate(self.strRefFixedRate,  xy=(0.70, 1.02), xycoords='axes fraction',size=10)
            else:
                ax.annotate(self.strTestFixedRate, xy=(0.70, 1.02), xycoords='axes fraction',size=10)


                
            bx = plt.gca()
            bx.format_coord = lambda x,y: '%12f, %12f' % (x,y)

            ax = plt.subplot(3,1,2)
            plt.plot(accutime,wholepitch, marker='.', linestyle='None', color=colorStyle, markersize=3,label=legendLabel, rasterized=True)
            plt.ylabel('[Deg]')
            plt.title('Pitch')
            plt.legend()
            bx = plt.gca()
            bx.format_coord = lambda x,y: '%12f, %12f' % (x,y)
            if bReference==False and self.bTestHeadingAvailable == False:
                if self.bTestGPSavailable:
                    plt.xlabel('Seconds of Week')
                else:
                    plt.xlabel('Seconds')

            if bReference==True or self.bTestHeadingAvailable == True:
                ax = plt.subplot(3,1,3)
                plt.plot(accutime,wholeyaw, marker='.', linestyle='None', color=colorStyle, markersize=3,label=legendLabel, rasterized=True)
                plt.ylabel('[Deg]')
                plt.title('Azimuth')
                plt.legend()
                bx = plt.gca()
                bx.format_coord = lambda x,y: '%12f, %12f' % (x,y)
                
                if bReference or self.bTestGPSavailable:
                    plt.xlabel('Seconds of Week')
                else:
                    plt.xlabel('Seconds')
            
            ax.annotate(self.timeStr, xy=(-0.1, -0.4), xycoords='axes fraction',size=9)
            if self.bGeneratePDF == True:
                pdf.savefig(fig)


            plt.close()

        return

    def normalangle(self,angle):
        ''' angle normalization
        '''
        a = angle
        while a<-180.0:
            a = a + 360.0
        while a>=180.0:
            a = a - 360.0
        return a

    def subangle(self,ang1,ang2):
        ''' angles substraction
        '''
        ang = ang1 - ang2
        ang = self.normalangle(ang)
        return ang

    def addangle(self,ang1,ang2):
        ''' angles addition
        '''
        ang = ang1 + ang2
        ang = self.normalangle(ang)
        return ang

    def getweeknum(self, weekseconds):
        ''' GPS week number from total GPS seconds
        '''
        return math.floor(weekseconds/(7*24*3600))

    def getweeksec(self, weekseconds):
        ''' return second of week from total GPS seconds
        '''
        return weekseconds - self.getweeknum(weekseconds)*(7*24*3600)

    def yearfour(self, year):
        ''' change two digit year to four digit
        '''
        if year<=80:
            year += 2000
        elif year<1990 and year>80:
            year += 1900
        return year

    def isleapyear(self, year):
        ''' find whether one year is a leap year
        '''
        return (self.yearfour(year)%4==0 and self.yearfour(year)%100!=0) or self.yearfour(year)%400==0
                     
    def timefromGPS(self, weeknum,weeksec):
        ''' return calendar time from GPS week number and SOW
        '''
        year = 0
        month = 0
        day = 0
        hour = 0
        minute = 0
        second = 0
        doy = 0
        daypermon = [31,28,31,30,31,30,31,31,30,31,30,31]

        weeknum += self.getweeknum(weeksec)
        weeksec  = self.getweeksec(weeksec)
        
        weekmin  = math.floor(weeksec/60.0)
        second   = weeksec - weekmin*60.0
        weekhour = math.floor(weekmin/60)
        minute   = weekmin - weekhour*60
        weekday  = math.floor(weekhour/24)
        hour     = weekhour - weekday*24

        totalday = weekday+weeknum*7
        if totalday<360:
            year = 1980
        else:
            year = 1981
            totalday -= 360
            while True:
                if totalday<365:
                    break
                if self.isleapyear(year): totalday -= 1
                totalday -= 365
                year += 1
        doy = totalday

        if totalday <= daypermon[0]:
            month = 1
        else:
            totalday -= daypermon[0]
            if self.isleapyear(year): totalday -= 1
            month = 2
            while True:
                if totalday<=daypermon[month-1]:
                    break
                else:
                    totalday -= daypermon[month-1]
                    month += 1
        if month==2 and self.isleapyear(year): totalday += 1
        day = math.floor(totalday)
        return [year,month,day,hour,minute,second,doy]

    def dayofyear(self, year, mon, day):
        ''' return day of year from calendar time
        '''
        daypermon = [31,28,31,30,31,30,31,31,30,31,30,31]
        totalDay  = day
        monIndex  = 0
        while monIndex<mon-1:
            totalDay += daypermon[monIndex]
            monIndex += 1

        if mon>2 and self.isleapyear(year):
            totalDay += 1
        return totalDay


    def convert2GPStime(self, year, month, day, hour, minute, second):
        ''' convert UTC calenday time to GPS Week, SOW, and DOY
            leap seconds have been corrected
        '''        
        if year<1981:
            totalDay = 0
        else:
            totalDay = 360

        yearIndex    = 1981
        while yearIndex<year:
            totalDay += 365
            if self.isleapyear(yearIndex):
                totalDay += 1
            yearIndex +=1

        doy = self.dayofyear(year, month, day)
        totalDay += doy
        weeknum  = math.floor(totalDay/7)
        weeksec  = (totalDay-weeknum*7)*24.0*3600+hour*3600+minute*60+second+self.leapseconds
        return [weeknum, weeksec, doy]

 #####       

if __name__ == "__main__":
    namestr1 = sys.argv[1]
    namestr2 = sys.argv[2]
    sysCompare = sysEvaluation()
    
    sysCompare.ReferenceTimeAlign      = float(sys.argv[3])
    sysCompare.TestTimeAlign           = float(sys.argv[4])
    sysCompare.initialTimeForRollPitch = float(sys.argv[5])
    sysCompare.initialTimeForYaw       = float(sys.argv[5])
    
    sysCompare.Open_Ref_Test_Files(namestr1)
    sysCompare.Open_Ref_Test_Files(namestr2)
    sysCompare.showmedifferences()
    sysCompare.showmeAttitudes()
    sysCompare.showmeAttitudesDifferences()




    
