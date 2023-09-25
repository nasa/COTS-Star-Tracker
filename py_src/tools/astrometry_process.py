#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
astrometry_process.py

This script is intended to process images through astrometry.net
in order to provide an independent attitude measurement
'''

################################
#LOAD LIBRARIES
################################
import os
import csv
import time
import pandas as pd
from datetime import datetime
from astroquery.astrometry_net import AstrometryNet

ast = AstrometryNet()
print("\n\n    !!!  Be sure to input your astrometry.net API key in the astrometry_process.py file  !!!!\n\n")

################################
#USER INPUT
################################
ast.api_key = 'YOUR_KEY_GOES_HERE' #astrometry.net API key
target_dir = "your/full/file/path/goes/here" #full path to the directory where the images to be processed are
solver_timeout = 3600 #seconds
crash_recovery = False #set this flag to True if the tool failed mid-run to resume where you left off

################################
#MAIN CODE
################################
#get list of images
the_files = os.listdir(target_dir)
image_list = []
for file in the_files:
    if (".png" in file) or (".jpg" in file) or (".jpeg" in file) or (".tiff" in file):
        image_list += [os.path.join(target_dir,file)]

iterator_image_list = image_list.copy() #required to unlink

#define vars
n = 0
ra = []
dec = []
solve_time = []
processed_image_list = []

#crash recovery
if crash_recovery:
    print("resuming from previous failed run")
    current_dir = os.getcwd()
    files = os.listdir(current_dir)
    biggest_file = ''
    biggest_file_size = 0
    for file in files:
        if 'astrometry_output.csv' in file:
            if os.path.getsize(os.path.join(current_dir,file)) > biggest_file_size:
                biggest_file = os.path.join(current_dir, file)
                biggest_file_size = os.path.getsize(os.path.join(current_dir,file)) 
    existing_data = pd.read_csv(biggest_file)
    already_solved_images = existing_data['image_name'].to_list()
    for image in already_solved_images:
        iterator_image_list.remove(image)
        n+=1

    processed_image_list = existing_data['image_name'].to_list()
    ra = existing_data['ra'].to_list()
    dec = existing_data['dec'].to_list()


for image_name in iterator_image_list:
    n+=1
    wcs_header = False
    try_again = True
    submission_id = None
    start_time = time.time()

    while try_again:
        try:
            if not submission_id:
                wcs_header = ast.solve_from_image(image_name, submission_id=submission_id, solve_timeout=solver_timeout, crpix_center=True)
            else:
                wcs_header = ast.monitor_submission(submission_id, solve_timeout=solver_timeout)
        except:# TimeoutError as e:
            #submission_id = e.args[1]
            print("\n Failure!  Maybe a timeout?")
            try_again = False
        else:
            # got a result, so terminate
            try_again = False

    print("\n"+str((n/len(image_list)*100))[:8]+"% complete overall")
    print(image_name)
    processed_image_list+=[image_name]
    if wcs_header:
        # Code to execute when solve succeeds
        solve_time += [time.time()-start_time]
        print("Found solution, took "+str(solve_time[-1])+" seconds\n\n")
        ra += [wcs_header['CRVAL1']]
        dec += [wcs_header['CRVAL2']]
    else:
        # Code to execute when solve fails
        solve_time += [time.time()-start_time]
        print("Unable to get solution, took "+str(solve_time[-1])+" seconds\n\n")
        ra += [999]
        dec += [999]

    #assemble dict
    the_data = {"image_name":processed_image_list,"ra":ra,"dec":dec,"solve_time":solve_time}
    now = str(datetime.now())
    now = now.split('.')
    now = now[0]
    now = now.replace(' ','_')
    now = now.replace(':','-')

    #write stuff
    keys=sorted(the_data.keys())
    with open(os.path.join(os.getcwd(),  now+'_astrometry_output.csv'),'w', newline='') as csv_file:
        writer=csv.writer(csv_file)
        writer.writerow(keys)
        writer.writerows(zip(*[the_data[key] for  key in keys]))

#assemble dict
the_data = {"image_name":processed_image_list,"ra":ra,"dec":dec,"solve_time":solve_time}
now = str(datetime.now())
now = now.split('.')
now = now[0]
now = now.replace(' ','_')
now = now.replace(':','-')

#write stuff
keys=sorted(the_data.keys())
with open(os.path.join(os.getcwd(),  now+'_astrometry_output_FINAL.csv'),'w', newline='') as csv_file:
    writer=csv.writer(csv_file)
    writer.writerow(keys)
    writer.writerows(zip(*[the_data[key] for  key in keys]))


print("\n\n      COMPLETE!\n\n")


