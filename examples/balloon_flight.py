#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""

This script is intended to capture imagery
and data from power-on throughout a
balloon flight on a raspberry pi

"""

################################
#LOAD LIBRARIES
################################
import os
import cv2
import csv
import json
import time
import psutil
import subprocess
import numpy as np
from datetime import datetime
from star_tracker import main
from star_tracker.cam_matrix import *
from star_tracker.array_transformations import *

################################
#USER INPUT
################################
starMatchPixelTol = 9
min_star_area = 3
max_star_area = 30
nmatch = 5
low_thresh_pxl_intensity = 15
hi_thresh_pxl_intensity = 255

VERBOSE = True
data_path = "data_path_goes_here"
cam_config_file = "xic_ximea_cam_example.json"

################################
#SUPPORT FUNCTIONS
################################
def init_ximea():
    #create instance for first connected camera 
    cam = xiapi.Camera()

    #start communication
    cam.open_device()

    #settings
    cam.set_exposure(50000) #microseconds
    cam.set_gain(24.0) #dB, max val varies per camera
    print('\nSet exposure to %s us \n' %cam.get_exposure())
    #looks like there are lots of auto-gain settings.  Unclear what to do with those.

    #create instance of Image to store image data and metadata
    img = xiapi.Image()

    #start data acquisition
    print('Starting data acquisition...')
    cam.start_acquisition()

    return cam, img




################################
#MAIN CODE
################################
#create directories
data_dir = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), 'data')

init_time = str(datetime.now())
init_time = init_time .split('.')
init_time = init_time [0]
init_time = init_time .replace(' ','_')
init_time = init_time .replace(':','-')
the_dir = os.path.join(data_dir,  init_time )
os.mkdir(the_dir)


if VERBOSE:
    print('\n data dir: ' + the_dir + '\n')

#save git state
home = os.getcwd()
os.chdir(os.path.dirname(os.path.abspath(__file__)))
git_show = str(subprocess.check_output(['git show'],shell=True))
git_status = str(subprocess.check_output(['git status'],shell=True))
git_dict = {'time':init_time,'git show':git_show,'git status':git_status}

with open(os.path.join(the_dir,'git_status.json'),"w") as statfile:
    json.dump(git_dict, statfile)

os.chdir(home)


#load star tracker stuff
cam_cal_dir = os.path.join(data_dir, 'cam_config')
input_data_dir = os.path.join(data_dir, data_path)
darkframe_file = 'xic_darkframe.jpg'
darkframe_file = os.path.join(input_data_dir, input_data_dir + darkframe_file)
print("darkframe: " + darkframe_file)
cam_file = os.path.join(cam_cal_dir, cam_config_file)


k = np.load(os.path.join(data_dir, 'k.npy'))
m = np.load(os.path.join(data_dir, 'm.npy'))
q = np.load(os.path.join(data_dir, 'q.npy'))
x_cat = np.load(os.path.join(data_dir, 'u.npy'))
indexed_star_pairs = np.load(os.path.join(data_dir, 'indexed_star_pairs.npy'))

camera_matrix, _, _ = read_cam_json(cam_file)
dx = camera_matrix[0, 0]
isa_thresh = starMatchPixelTol*(1/dx)


#adjust USB buffer size (required for systems like Pis)
if VERBOSE:
    os.system('cat /sys/module/usbcore/parameters/usbfs_memory_mb')
os.system("/bin/bash -c 'tee /sys/module/usbcore/parameters/usbfs_memory_mb >/dev/null <<<0'")
if VERBOSE:
    os.system('cat /sys/module/usbcore/parameters/usbfs_memory_mb')


#initialize cameras
from ximea import xiapi
cam, img = init_ximea()


#create variables for data logging
time_interval =   20 #sec

ttime = []
stemp = []
sram  = []
scpu  = []
image_time = []
solve_time = []
qs = []
qv0 = []
qv1 = []
qv2 = []
start_time = time.time()

#capture imagery and data forever
n = 0
while True:

    #reset time
    end_time = time.time()

    #get data and pass them from camera to img

    try:

        n+=1
        if VERBOSE:
            print('capture ' + str(n))

        gain_list = [0.0,12.0,24.0]#24.0

        for gain in gain_list:

            cam.set_exposure(50000) #microseconds
            cam.set_gain(gain) #dB, max val varies per camera
            print('\nSet exposure to %s us \n' %cam.get_exposure())
            time.sleep(1)
            img_start_time = time.time()
            cam.get_image(img)

            image_time += [time.time()-img_start_time]

            #create numpy array with data from camera. Dimensions of array are determined
            #by imgdataformat
            frame = img.get_image_data_numpy()

            #save image
            now = str(datetime.now())
            now = now.split('.')
            now = now[0]
            now = now.replace(' ','_')
            now = now.replace(':','-')
            image_filename = os.path.join(data_dir, init_time, str(n)+'_'+now+"_50ms_"+str(gain)+".jpg")
            cv2.imwrite(image_filename,frame)

            #save images at varying exposure
            cam.set_exposure(100000) #microseconds
            cam.set_gain(gain) #dB, max val varies per camera
            print('Set exposure to %s us' %cam.get_exposure())
            time.sleep(1)
            cam.get_image(img)
            frame2 = img.get_image_data_numpy()
            image_filename2 = os.path.join(data_dir, init_time, str(n)+'_'+now+"_100ms_"+str(gain)+".jpg")
            cv2.imwrite(image_filename2,frame2)

            cam.set_exposure(150000) #microseconds
            cam.set_gain(gain) #dB, max val varies per camera
            print('Set exposure to %s us' %cam.get_exposure())
            time.sleep(1)
            cam.get_image(img)
            frame2 = img.get_image_data_numpy()
            image_filename2 = os.path.join(data_dir, init_time, str(n)+'_'+now+"_150ms_"+str(gain)+".jpg")
            cv2.imwrite(image_filename2,frame2)

            cam.set_exposure(200000) #microseconds
            cam.set_gain(gain) #dB, max val varies per camera
            print('Set exposure to %s us' %cam.get_exposure())
            time.sleep(1)
            cam.get_image(img)
            frame2 = img.get_image_data_numpy()
            image_filename2 = os.path.join(data_dir, init_time, str(n)+'_'+now+"_200ms_"+str(gain)+".jpg")
            cv2.imwrite(image_filename2,frame2)

            cam.set_exposure(250000) #microseconds
            cam.set_gain(gain) #dB, max val varies per camera
            print('Set exposure to %s us' %cam.get_exposure())
            time.sleep(1)
            cam.get_image(img)
            frame2 = img.get_image_data_numpy()
            image_filename2 = os.path.join(data_dir, init_time, str(n)+'_'+now+"_250ms_"+str(gain)+".jpg")
            cv2.imwrite(image_filename2,frame2)


    except:
        faill = True
        while faill == True:
            print('\n\n\nFAILURE!  Attempting to recover...\n\n\n')
            time.sleep(2)
            try:
                cam.stop_acquisition()
                cam.close_device()
                print('device closed')
            except:
                print('device close failed')
                pass

            try:
                cam, img = init_ximea()
                cam.get_image(img)
                faill = False
                n+=1
            except:
                print('FDIR re-init failed')
                pass

        #create numpy array with data from camera. Dimensions of array are determined
        #by imgdataformat
        frame = img.get_image_data_numpy()

        #save image
        now = str(datetime.now())
        now = now.split('.')
        now = now[0]
        now = now.replace(' ','_')
        now = now.replace(':','-')
        image_filename = os.path.join(data_dir, init_time, str(n)+'_'+now+"_50ms.jpg")
        cv2.imwrite(image_filename,frame)



    #run star tracker
    solve_start_time = time.time()

    q_est, idmatch, nmatches, x_obs, rtrnd_img = main.star_tracker(
            image_filename, cam_file, darkframe_file=darkframe_file, m=m, q=q, x_cat=x_cat, k=k, indexed_star_pairs=indexed_star_pairs, graphics=False,
            min_star_area=min_star_area, max_star_area=max_star_area, isa_thresh=isa_thresh, nmatch=nmatch)

    solve_time += [time.time()-solve_start_time]

    #collect data
    try:
        assert not np.any(np.isnan(q_est))
        if VERBOSE:
            print('est q: ' + str(q_est))
        qs += [q_est[3]]
        qv0 += [q_est[0]]
        qv1 += [q_est[1]]
        qv2 += [q_est[2]]
    except AssertionError:
        if VERBOSE:
            print('NO VALID STARS FOUND\n')
        qs += [999]
        qv0 += [999]
        qv1 += [999]
        qv2 += [999]


    temp = subprocess.check_output(['sudo', 'cat', '/sys/class/thermal/thermal_zone0/temp'])

    ttime += [time.time()]
    stemp += [float(temp)/1000]
    sram  += [psutil.virtual_memory().percent]
    #scpu  += [psutil.cpu_percent(2)]
    scpu  += [psutil.cpu_percent()]

    if end_time-start_time > time_interval:

        data = {'time':ttime,'system temp (C)':stemp,'RAM':sram,'CPU':scpu, 'image capture time (s)': image_time,
                'image solve time (s)':solve_time, 'qs':qs,'qv0':qv0,'qv1':qv1,'qv2':qv2}

        now = str(datetime.now())
        now = now.split('.')
        now = now[0]
        now = now.replace(' ','_')
        now = now.replace(':','-')

        #write stuff
        keys=sorted(data.keys())
        with open(os.path.join(data_dir,  init_time , now+'sys_data_file.csv'),'w') as csv_file:
             writer=csv.writer(csv_file)
             writer.writerow(keys)
             writer.writerows(zip(*[data[key] for  key in keys]))

        #reset vars
        ttime = []
        stemp = []
        sram  = []
        scpu  = []
        image_time = []
        solve_time = []
        qs = []
        qv0 = []
        qv1 = []
        qv2 = []

        #force changes out of memory and onto disc
        if VERBOSE:
            print('sync')
        os.system('sync')

        #reset time
        start_time = time.time()




#stop data acquisition
print('Stopping acquisition...')
cam.stop_acquisition()

#stop communication
cam.close_device()

# Destroys the OpenCv windows
cv2.destroyAllWindows()


print("THE END")





