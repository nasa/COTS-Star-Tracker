#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
pi_cam_interface.py

This script is intended to provide a 
simple interface to a raspberry pi camera

"""

################################
#LOAD LIBRARIES
################################
import os
from picamera import PiCamera

################################
#USER INPUT
################################
#select camera type
pi_cam = "v2" # v1, v2, hq
image_extension = ".jpg"

################################
#MAIN CODE
################################
#create camera object
camera = PiCamera()

#camera options
v1_res = [2592, 1944]
v2_res = [3280, 2464]

if pi_cam == "v1":
    selected_cam_res = v1_res
if pi_cam == "v2":
    selected_cam_res = v2_res
if pi_cam == "hq":
    selected_cam_res = v2_res #TODO update to actual HQ res, but somehow fix memory issues

#set cam parameters
camera.resolution = (selected_cam_res)
#camera.brightness = 50  #does not appear useful
#camera.contrast = 50 #may want this high
#camera.awb_mode = 'sunlight'  #this sets color balance
#camera.exposure_mode = 'night' #this drives gain settings
camera.iso=800 #set max iso

camera.start_preview(alpha=200, resolution=(int(selected_cam_res[0]/3),int(selected_cam_res[1]/3)))

usr_in = ''
capture_number = 0

while usr_in == '':
    
    usr_in = input('Press enter to capture, input any character to exit.\n')
    
    capture_name = 'pi_test_'+str(capture_number)+image_extension
    
    #Rename until name is not taken
    while os.path.exists(os.path.join(os.getcwd(), capture_name)):
        capture_number += 1
        capture_name = 'pi_test_'+str(capture_number)+'.jpg'
                
    if usr_in == '':
        camera.capture(os.path.join(os.getcwd(),capture_name))
        print(capture_name+' saved!\n')

#release cam resources, prevents GPU memory leaks
camera.close()

