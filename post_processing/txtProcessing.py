import matplotlib.pyplot as plt
import numpy as np
import sys
import csv
import requests
import time
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
from matplotlib.backends.backend_pdf import PdfPages

from azure.storage.blob import AppendBlobService
from azure.storage.blob import ContentSettings
from azure.storage.blob import BlockBlobService

raw_data = []

time = []
roll = []
pitch = []
wx = []
wy = []
wz = []
lat = []
long = []
cog = []
sog = []
yaw_gps = []
pitch_gps = []
roll_gps = []
ax = []
ay = []
az = []

for line in urlopen('https://navview.blob.core.windows.net/data-test/'+ sys.argv[1] +'.txt'):
    line = line.strip()
    raw_data.append(line.decode('utf-8').split(','))

length = len(raw_data)  
for i in range(1,length):
    time.append(raw_data[i][0])
    roll.append (raw_data[i][1])
    pitch.append(raw_data[i][2])
    wx.append (raw_data[i][3])
    wy.append(raw_data[i][4])
    wz.append (raw_data[i][5])
    lat.append(raw_data[i][6])
    long.append (raw_data[i][7])
    cog.append(raw_data[i][8])
    sog.append (raw_data[i][9])
    yaw_gps.append(raw_data[i][10])
    pitch_gps.append (raw_data[i][11])
    roll_gps.append(raw_data[i][12])
    ax.append (raw_data[i][13])
    ay.append(raw_data[i][14])
    # az.append (raw_data[i][15])
 
for i in range(2,len(time)):
    # print()
    time[i] = float(time[i]) - float(time[1])
    
time[1] = 0.0     
print(time)
data = [time,roll,pitch,wx,wy,wz,lat,long,cog,sog,yaw_gps,pitch_gps,roll_gps,ax,ay]

# for i in range(1,len(data)):
#         if raw_data[i] != 'time':
#             plt.figure(figsize=(9, 9))
#             plt.plot(range(len(time)),data[i], label='Loaded from file!')
#             plt.xlim([0,len(time)]) # set x value range
#             # plt.xticks(np.arange(0, len(time), step=len(time)/10)) 
#             plt.xticks(np.arange(0, len(time), 10))
#             # plt.xticklabels('auto')
#             # plt.set_xtickslabels(x[::10])
#             plt.title(raw_data[0][i])
            # plt.savefig(str(raw_data[0][i]) + ".png")

with PdfPages('./post_processing/data/'+ sys.argv[1] +'_plots.pdf') as pdf:
    for i in range(1,len(data)):
        if raw_data[i] != 'time':
            plt.figure(figsize=(9, 9))
            plt.plot(range(len(time)),data[i], label='Loaded from file!')
            plt.xlim([0,len(time)]) # set x value range
            # plt.xticks(np.arange(0, len(time), step=len(time)/10)) 
            plt.xticks(np.arange(0, len(time), 10))
            # plt.xticklabels('auto')
            # plt.set_xtickslabels(x[::10])
            plt.title(raw_data[0][i])
            # plt.savefig(str(raw_data[0][i]) + ".png")
            pdf.savefig()  # saves the current figure into a pdf page
    # plt.show()      
    plt.close()

local_file_name = './post_processing/data/'+ sys.argv[1] +'_plots.pdf'
# Create a blob client using the local file name as the name for the blob
block_blob_service = BlockBlobService(account_name='navview',
                                      account_key='hiIfcpU8zh75BbbOWEY5cZLyjc/DZZq4zoJuLUToSd3Wphi0Lc4upSOg3ZDh0RTZ5lhnrXXhTRoONPl2HRfrXQ==', 
                                    #   account_key=os.getenv("ACCOUNT_KEY"), # account_key to be encripted
                                      protocol='http')
block_blob_service.create_blob_from_path('data-test', sys.argv[1] +'_plots.pdf', local_file_name)

time.sleep(1)
response = requests.get("http://localhost:3060/api/drive/complete/")



    