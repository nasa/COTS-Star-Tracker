#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
create_darkframe.py

Script used for creating dark frame from existing set of
star images

"""

################################
#LOAD LIBRARIES
################################
import os
import cv2
import time
import star_tracker.ground as ground


#############################################################
#USER-INPUT VARS
#############################################################
path_to_images = "path/to/your/images" #full path to the directory where the images to be processed are
num_images = 7 #default of 7


#############################################################
#MAIN CODE
#############################################################

# CREATE DARK FRAME
print("\ncreating dark frame...")
start_time = time.time()

#img_list = ground.find_files_pattern(im_pattern, im_path, exclude='dark')
img_list = [os.path.join(path_to_images,f) for f in os.listdir(path_to_images) if os.path.isfile(os.path.join(path_to_images,f))]

if len(img_list) < 1:
    print("No images found, unable to create darkframe")

else:
    darkframe = ground.create_darkframe(img_list, num_images)

    # Save dark frame
    darkframe_path = os.path.join(path_to_images,'autogen_darkframe.jpg')
    cv2.imwrite(darkframe_path, darkframe)
    print(darkframe_path)
    print("\n...dark frame creation complete in {0} seconds\n".format(time.time()-start_time))


