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
import time
import star_tracker.ground as ground
import cv2


#############################################################
#USER-INPUT VARS
#############################################################
#data_path = "data/StarFields/"
#im_path = "Stars_250ms/"

#data_path = "data/synth_stars/"
#im_path = "Pass_F_WOut-Dist_WOut-SA/"

#data_path = "data/StarsNY/"
#im_path = "Up_North/"


#############################################################
#MAIN CODE
#############################################################

tools_dir = os.path.dirname(os.path.realpath(__file__))
py_src_dir = os.path.dirname(tools_dir)
repo_dir = os.path.dirname(py_src_dir)
cam_config_dir = os.path.join(repo_dir, os.path.join('data','cam_config'))


# CREATE DARK FRAME
print("\ncreating dark frame...")
start_time = time.time()

num_images = 7
#img_list = ground.find_files_pattern(im_pattern, im_path, exclude='dark')
darkframe_dir = os.path.join(repo_dir, data_path)
df_config = os.path.join(darkframe_dir, im_path)
img_list = [os.path.join(df_config,f) for f in os.listdir(df_config) if os.path.isfile(os.path.join(df_config,f))]

if len(img_list) < 1:
    print("No images found, unable to create darkframe")

else:
    darkframe = ground.create_darkframe(img_list, num_images)

    # Save dark frame
    darkframe_path = os.path.join(df_config,'autogen_darkframe.jpg')
    cv2.imwrite(darkframe_path, darkframe)
    print(darkframe_path)
    print("\n...dark frame creation complete in {0} seconds\n".format(time.time()-start_time))


