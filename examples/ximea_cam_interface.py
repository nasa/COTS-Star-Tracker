#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""

This script is intended to provide a 
quick and dirty interface to a
ximea camera

"""

################################
#LOAD LIBRARIES
################################
import os
import cv2
import time
import numpy as np
from ximea import xiapi


################################
#USER INPUT
################################
exposure = 20000 #microseconds
gain = 24.0 #scale varies by camera (0-18 or 0-24)
image_extension = ".jpg"


################################
#MAIN CODE
################################
#adjust buffer size
os.system('cat /sys/module/usbcore/parameters/usbfs_memory_mb')
os.system("/bin/bash -c 'tee /sys/module/usbcore/parameters/usbfs_memory_mb >/dev/null <<<0'")
os.system('cat /sys/module/usbcore/parameters/usbfs_memory_mb')


#create instance for first connected camera 
cam = xiapi.Camera()

#start communication
print('Opening first camera...')
cam.open_device()

#settings
cam.set_exposure(exposure)
cam.set_gain(gain) 
#looks like there are lots of auto-gain settings.  Unclear what to do with those.

#create instance of Image to store image data and metadata
img = xiapi.Image()

#start data acquisition
print('Starting acquisition...')
cam.start_acquisition()

#display image
n = 0
while True:

    #get data and pass them from camera to img
    try:
        cam.get_image(img)
    except:
        time.sleep(1)
        print('failure-- resetting device connection')
        cam.stop_acquisition()
        cam.close_device()
        try:
            cam = xiapi.Camera()
            cam.open_device()
            cam.set_exposure(20000)
            cam.set_gain(24.0) #scale from 0-18??
            img = xiapi.Image()
            cam.start_acquisition()

            time.sleep(1)
            cam.get_image(img)
        except:
            pass

    #create numpy array with data from camera. Dimensions of array are determined
    #by imgdataformat
    array = img.get_image_data_numpy()

    # ...reshape it in an numpy array...
    frame = array

    # ...resize the image by a half
    frame2 = cv2.resize(frame,(0,0),fx=0.5, fy=0.5)

    #...and finally display it
    cv2.imshow("SimpleLive_Python_XIMEA_OpenCV", frame2)

    # Press q if you want to end the loop and p to take a picture
    if cv2.waitKey(1) & 0xFF == ord('p'):
        cv2.imwrite("ximea_"+str(n)+image_extension,frame)
        n+=1
    elif cv2.waitKey(1) & 0xFF == ord('q'):
        cv2.imwrite("ximea_"+str(n)+image_extension,frame)
        break


#stop data acquisition
print('Stopping acquisition...')
cam.stop_acquisition()

#stop communication
cam.close_device()

# Destroys the OpenCv windows
cv2.destroyAllWindows()


print("THE END")
