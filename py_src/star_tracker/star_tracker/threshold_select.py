#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
threshold_select.py

This program is intended to:
* open an image
* contour bright spots
* centroid contours
* undistort centroids
* output corresponding unit vectors in the camera frame

'''


################################
#LOAD LIBRARIES
################################
import os
import sys
import cv2
import json
import time
#from math import *
#from PIL import Image
import matplotlib as mpl
mpl.use('tkagg')  #set this backend for rpi
import matplotlib.pyplot as plt
import numpy as np
from cam_matrix import read_cam_json

################################
#USER INPUT
################################
#define inputs
#im_filename = "../data/IDS-HWIL/ids_0.jpg"
#cam_config = "IDS_cam_4real.json"
#darkframe_file = os.path.join('../data/IDS-HWIL/IDS_darkframe.jpg')
im_filename = "../data/ximea-HWIL/ximea_8.jpg"
cam_config = "ximea_cam.json"
darkframe_file = os.path.join('../data/ximea-HWIL/ximea_darkframe_20ms_18g.jpg')

thresh_hi = 255
thresh_lo = 60 #127

cam_config_dir = os.path.join(os.path.dirname(os.getcwd()), 'data','cam_config')
cam_config_file = os.path.join(cam_config_dir, cam_config)


################################
#MAIN CODE
################################
if os.path.exists(cam_config_file):
    print("reading cam params from file")
    camera_matrix, cam_resolution, dist_coefs = read_cam_json(cam_config_file)

else:
    print(cam_config_file)
    print("file not found, manually defining cam params")

    cam_resolution = [3280, 2464] #pixels
    pixel_pitch = [0.00112, 0.00112] #mm
    focal_length = 3 #mm
    prncple_pt = [cam_resolution[0]/2,cam_resolution[1]/2] #pixels
    k_1 = 0
    k_2 = 0
    p_1 = 0
    p_2 = 0
    k_3 = 0
    k_4 = 0
    k_5 = 0
    k_6 = 0

    #convert to opencv units
    prncple_pt_cv2 = [prncple_pt[0]*pixel_pitch[0],prncple_pt[1]*pixel_pitch[1]] #mm
    focal_length_cv2 = [focal_length/pixel_pitch[0],focal_length/pixel_pitch[1]] #pixels

    #populate camera matrix
    camera_matrix = np.array([[focal_length_cv2[0], 0,        prncple_pt_cv2[0]],
                         [0,         focal_length_cv2[1], prncple_pt_cv2[1]],
                         [0,                  0,          1]])

    #distortion coeffs
    dist_coefs =  np.array([k_1, k_2, p_1, p_2, k_3, k_4, k_5, k_6])

#start the clock and load an image
start_time = time.time()
start_time1 = time.time()
img = cv2.imread(im_filename)
#img = cv2.imread(im_filename, cv2.IMREAD_GRAYSCALE)  #this causes purple images???
print("\nImage read in " + str(time.time()-start_time) + "s ")
plt.imshow(img,cmap='Greys_r'), plt.show()
start_time = time.time()


#load darkframe (if it exists), convert, subtract
try:
    darkframe = cv2.imread(darkframe_file)
    img2 = cv2.subtract(img, darkframe)
    img = img2
except:
    print('\n\ndark frame failed!\n\n')


# convert the image to grayscale
gray_image = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
#gray_image = img #TROUBLESHOOTING
#plt.imshow(gray_image,cmap='Greys_r'), plt.show()

# convert the grayscale image to binary image
ret, thresh = cv2.threshold(gray_image, thresh_lo, thresh_hi, cv2.THRESH_BINARY)


# find contours in the binary image
try:
    contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
except:
    junk_im, contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
print("\nContours detected in " + str(time.time()-start_time) + "s ")

print('darkframe')
plt.imshow(darkframe), plt.show()
print('subtracted image')
plt.imshow(img2), plt.show()
print('thresholded image')
plt.imshow(thresh,cmap='Greys_r'), plt.show()



start_time = time.time()

#centroid blobs
found_centroids = []
print("\nfound " + str(len(contours)) + " contours")
for c in contours:
   # calculate moments for each contour
   M = cv2.moments(c)

   # calculate x,y pixel location of center
   if M["m00"] > 0:
      cX = M["m10"] / M["m00"]
      cY = M["m01"] / M["m00"]

      found_centroids += [[cX,cY]]

   else: #if there's not enough of an area to bound, use sole point?  Look more into this...
      cX = c[0][0][0]
      cY = c[0][0][1]

      found_centroids += [[float(cX),float(cY)]] #as it turns out, undistort hates ints

   cv2.circle(img, (int(cX), int(cY)), 2, (0, 255, 0), -1)


print("\nCentroids created in " + str(time.time()-start_time) + "s ")
start_time = time.time()


print("\nfound " + str(len(found_centroids)) + " centroids")
if len(found_centroids) < 1:
    print('\n\n\nNO CENTROIDS FOUND.  Quitting.\n\n')
    sys.exit()
found_centroids = np.array([found_centroids])
found_centroids.reshape(len(found_centroids[0]),1,2)
#print('found centroids ' + str(found_centroids))




#undistort points
undistorted_centroids = cv2.undistortPoints(found_centroids, camera_matrix, dist_coefs)

print("\nCentroids undistorted in " + str(time.time()-start_time) + "s ")
start_time = time.time()

print("")
#print(undistorted_centroids)
print("")


#convert to unit vectors in camera frame
output_vecs = []

for centroid in undistorted_centroids[0]:

   the_array = np.array([centroid[0],centroid[1],1])
   vec_norm = np.linalg.norm([the_array])
   output_vecs += [[the_array[0]/vec_norm,the_array[1]/vec_norm,the_array[2]/vec_norm]]

print("\nvectors created in " + str(time.time()-start_time) + "s ")
start_time = time.time()

#print('\ncam vecs: ' + str(output_vecs))

print("\ntotal time " + str(time.time()-start_time1) + "s ")


# display the image
print('processed image w/ overlaid centroids')
plt.imshow(img), plt.show()

print('Done')


