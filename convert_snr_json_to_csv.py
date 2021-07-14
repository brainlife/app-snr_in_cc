#!/usr/bin/env python3

import os,sys
import json
import pandas as pd
import numpy as np

# load snr.json
with open('product.json','r') as snr_f:
	snr_data = json.load(snr_f)

# set up output dataframe
out_data = pd.DataFrame(columns={'snr','volumes','dir_x','dir_y','dir_z'})

# grab snr data in all directions
snr = []
snr = [ float(f.split(', ')[1]) for f in snr_data['SNR data in all directions'] ]
snr = snr + [ float(f) for f in snr_data['SNR in b0, X, Y, Z'][1:]]

# grab volumes data
volumes = []
volumes = [ f.split(', ')[0] for f in snr_data['SNR data in all directions'] ]
volumes = volumes + ['b0_X','b0_Y','b0_Z']

# grab the dir_x,y,z data
dir_x = []
dir_y = []
dir_z = []
dir_x = [ float(snr_data['b0, X, Y, Z directions'][1]) ]
dir_y = [ float(snr_data['b0, X, Y, Z directions'][2]) ]
dir_z = [ float(snr_data['b0, X, Y, Z directions'][3]) ]

dir_x = dir_x + [ float(f.split(', ')[1].replace('[ ','').replace(']','').replace('[','').split()[0]) for f in snr_data['direction vectors'] ]
dir_y = dir_y + [ float(f.split(', ')[1].replace('[ ','').replace(']','').replace('[','').split()[1]) for f in snr_data['direction vectors'] ]
dir_z = dir_z + [ float(f.split(', ')[1].replace('[ ','').replace(']','').replace('[','').split()[2]) for f in snr_data['direction vectors'] ]

dir_x = dir_x + [float(snr_data['b0, X, Y, Z directions'][1]),0,0]
dir_y = dir_y + [0,float(snr_data['b0, X, Y, Z directions'][2]),0]
dir_z = dir_z + [0,0,float(snr_data['b0, X, Y, Z directions'][3])]

# set columns in output dataframe
out_data['snr'] = snr
out_data['volumes'] = volumes
out_data['dir_x'] = dir_x
out_data['dir_y'] = dir_y
out_data['dir_z'] = dir_z

# reorder columns
out_data = out_data[['volumes','dir_x','dir_y','dir_z','snr']]

# save data
if ~os.path.isdir('./snr-stats'):
	os.mkdir('./snr-stats')

out_data.to_csv('./snr-stats/snr.csv',index=False)
