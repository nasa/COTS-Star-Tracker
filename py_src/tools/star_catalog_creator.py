#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
star_catalog_creator.py

Script used for creating star catalogs, k-vectors, etc.

Distributed under the 3-Clause BSD License (below)

Copyright 2019 Rensselaer Polytechnic Institute
(Dr. John Christian, Devin Renshaw, Grace Quintero)

Redistribution and use in source and binary forms,
with or without modification, are permitted provided
 that the following conditions are met:

1. Redistributions of source code must retain the above
 copyright notice, this list of conditions and the
 following disclaimer.

2. Redistributions in binary form must reproduce the above
 copyright notice, this list of conditions and the following
 disclaimer in the documentation and/or other materials
 provided with the distribution.

3. Neither the name of the copyright holder nor the names
 of its contributors may be used to endorse or promote
 products derived from this software without specific
 prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

################################
#LOAD LIBRARIES
################################
import os
import cv2
import time
import numpy as np
import star_tracker.ground as ground
import star_tracker.cam_matrix as cam_matrix

#############################################################
#USER-INPUT VARS
#############################################################
b_thresh = 6.0  # brightness threshold, recommend 5-6 to start with
cam_config_file = '' #the name of the camera config file in /data


#############################################################
#MAIN CODE
#############################################################
# CREATE CATALOG
#row_start = 50  # row where data begins in starCatParallax
save_vals = True  # option to save values or just return them from function
#col_brightness = 1  # column in star catalog containing the brightness values
#col_rade = [5, 6]   # column containing RA and DE values in starCatParallax
#col_pm = [3, 4]   # column containing proper motion values in starCatParallax
#col_par = 2   # column containing parallax values in starCatParallax
#row_ep = 12   # row with catalog epoch
#row_ep = 38   # row with catalog epoch


#Excess rows to remove from starcat_file
excess_rows = [53, 54]
# column (0-indexing) containing the Hipparcos ID number
index_col = 2


tools_dir = os.path.dirname(os.path.realpath(__file__))
py_src_dir = os.path.dirname(tools_dir)
repo_dir = os.path.dirname(py_src_dir)
starcat_file  = os.path.join(repo_dir, os.path.join('data','starcat.tsv'))
cam_config_dir = os.path.join(repo_dir, os.path.join('data','cam_config'))
cam_config_file = os.path.join(cam_config_dir, cam_config_file)

camera_matrix, cam_resolution, dist_coefs = cam_matrix.read_cam_json(cam_config_file)

nrow = cam_resolution[1]
ncol = cam_resolution[0]
fov = cam_matrix.cam2fov(cam_matrix.cam_matrix_inv(camera_matrix), nrow, ncol)

#TODO: Replace with actual observer location. Currently set to zero to just ignore parallax.
rB_scalar = 0*149597870.693
rB = np.array([[rB_scalar],[0],[0]]) #unit vector from SSB to observer

#check to see if the user has defined a directory to save things in already
try: print("Creating star pair catalog including stars up to mag "+str(b_thresh)+", saving to: "+str(save_dir)+" ...")
except: #if not, put them in the data dir
    file_dir = os.path.realpath(__file__)
    tools_dir = os.path.dirname(file_dir)
    py_src_dir = os.path.dirname(tools_dir)
    top_level_dir = os.path.dirname(py_src_dir)
    save_dir = os.path.join(top_level_dir,'data')
    print("Creating star pair catalog including stars up to mag "+str(b_thresh)+", saving to: "+str(save_dir)+" ...")
start_time = time.time()
ground.create_star_catalog(starcat_file=starcat_file, brightness_thresh=b_thresh,
                            excess_rows=excess_rows, index_col=index_col, fov=fov,
                            save_vals=save_vals, rB=rB, save_dir=save_dir)

print("\n...catalog creation complete in " + str(time.time()-start_time)+ " seconds\n")

