#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
checkerboard_cam_cal.py

This program is intended to:
* ingest images of checkerboard targets
* find the corners of the checkerboards
* generate camera cal parameters
* refine said parameters
* demo undistortion on the first loaded image
* save off camera params

'''

################################
#LOAD LIBRARIES
################################
import os
import sys
import cv2
import glob
import time
import json
import numpy as np

################################
#USER INPUT
################################
number_of_target_columns = 9 #9
number_of_target_rows = 6 #6
img_extension = '.jpg'

################################
#MAIN CODE
################################
# termination criteria
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

# prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
objp = np.zeros((number_of_target_rows*number_of_target_columns,3), np.float32)
objp[:,:2] = np.mgrid[0:number_of_target_columns,0:number_of_target_rows].T.reshape(-1,2)

# Arrays to store object points and image points from all the images.
objpoints = [] # 3d point in real world space
imgpoints = [] # 2d points in image plane.

images = glob.glob('*'+img_extension)

print("\n\nAssuming checkerboard is "+str(number_of_target_rows)+" by "+str(number_of_target_columns))
print('NOTE: ENSURE NUMBER OF SQUARES IS SET CORRECTLY (n-1 rows and n-1 columns)')

print("\nLooking for images ending in " + str(img_extension)+"...")
print('...found ' + str(len(images)) + ' images:')
print(images)

if len(images) == 0:
    print("""\nPlease ensure calibration images are present in the current
working directory and also have the following filetype(case sensitive):""")
    print("\n"+img_extension)
    print("\nExiting...")
    sys.exit()

elif len(images) < 10:
    print("""\nIt's recommended to use at least 10 images. Please 
ensure your lower image count is intentional, otherwise it is
recommended to add additional calibration images for improved results.""")
    usr_in = input("\nPress enter to continue.")


scale_percent = 100 # percent of original size

for fname in images:

    print("\n" + str(fname))

    img = cv2.imread(fname)

    width = int(img.shape[1] * scale_percent / 100)
    height = int(img.shape[0] * scale_percent / 100)
    dim = (width, height)
    img = cv2.resize(img, dim)

    gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)


    # Find the chess board corners
    print('Finding corners...')
    start_time = time.time()
    ret, corners = cv2.findChessboardCorners(gray, (number_of_target_columns,number_of_target_rows),None)
    print("...findChessboardCorners finished in " + str(time.time()-start_time) + "s")


    # If found, add object points, image points (after refining them)
    if ret == True:

        print("    found corners!")
        image_size = (img.shape[1],img.shape[0]) # tuple required, input args flipped
        objpoints.append(objp)

        cv2.cornerSubPix(gray,corners,(11,11),(-1,-1),criteria)
        imgpoints.append(corners)

        # Draw and display the corners
        cv2.drawChessboardCorners(img, (number_of_target_columns,number_of_target_rows), corners, ret)
        cv2.imshow('img',img)
        cv2.waitKey(500)


cv2.destroyAllWindows()


# THE CAL
start_time = time.time()
RMS_reproj_err_pix, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(objpoints, imgpoints, image_size, None, None, flags=(cv2.CALIB_FIX_PRINCIPAL_POINT))
[fovx, fovy, fl, pp, ar] = cv2.calibrationMatrixValues(mtx, image_size, 4.8e-3*image_size[0], 4.8e-3*image_size[1] ) # Note that this may not work without an accurate detector size (Args 3 and 4)

print("\nCamera parameter estimation (calibration) complete in " + str(time.time()-start_time) + "s")
print("\nReprojection Error (pix, RMS): " + str(RMS_reproj_err_pix))
#print("\nFoV (x,y): "+str(fovx)+", "+str(fovy))
#print("Focal Length (pix): " + str(flpix))
print("Focal Length (mm): " + str(fl))
print("\nCamera Matrix: "+ str(mtx))
print("\nDist: " + str(dist))
print("\nRvecs: " + str(rvecs))
print("\nTvecs: " + str(tvecs))
print("")
    

# undistort and remap to verify
print("\n\n"+str(images[0]) + " now being undistorted to demonstrate calibration.")
img = cv2.imread(images[0])

width = int(img.shape[1] * scale_percent / 100)
height = int(img.shape[0] * scale_percent / 100)
dim = (width, height)
img = cv2.resize(img, dim)

h,  w = img.shape[:2]
h2,  w2 = img.shape[:2]
newcameramtx, roi=cv2.getOptimalNewCameraMatrix(mtx,dist,(w,h),1,(w,h))

# undistort and save
dst = cv2.undistort(img, mtx, dist, None, newcameramtx)
cv2.imwrite('calibresult.png',dst)

# crop and save
x,y,w,h = roi
dst = dst[y:y+h, x:x+w]
cv2.imwrite('calibresult_cropped.png',dst)
print("...complete.  Please evaluate calibresult.png and calibresult_cropped.png to determine if they are less distorted than the original.")


# populate dict
dist_l=dist.tolist()
#newcameramtx_l = mtx.tolist()
newcameramtx_l = newcameramtx.tolist()
# in this case, cy is vp and cx is up.
cam_cal_dict = {'camera_matrix': newcameramtx_l, 'dist_coefs': dist_l, 'resolution':[image_size[0],image_size[1]], 'camera_model':'Brown','k1':dist_l[0][0],'k2':dist_l[0][1],'k3':dist_l[0][4],'p1':dist_l[0][2],'p2':dist_l[0][3],'fx':newcameramtx_l[0][0],'fy':newcameramtx_l[1][1],'up':newcameramtx_l[0][2],'vp':newcameramtx_l[1][2],'skew':newcameramtx_l[0][1], 'RMS_reproj_err_pix':RMS_reproj_err_pix}

#save if it's good
usr_in = input('\n\n\nIf you are satisfied with the results, provide a file name.  Otherwise, enter "no". [no/file_name]: ')

usr_in=usr_in.lower()
if usr_in in ['n', 'no']:
    sys.exit()
    
#save
if usr_in == '':
    usr_in = 'generic_cam_params'
usr_in_split = usr_in.split('.json')
usr_in = usr_in_split[0]
checkerboard_dir = os.path.dirname(os.path.realpath(__file__))
cam_cal_dir = os.path.dirname(checkerboard_dir)
tools_dir = os.path.dirname(cam_cal_dir)
py_src_dir = os.path.dirname(tools_dir)
repo_dir = os.path.dirname(py_src_dir)
cam_config_dir = os.path.join(repo_dir,'data','cam_config')
full_cam_file_path = os.path.join(cam_config_dir,usr_in+'.json')
with open(full_cam_file_path, 'w') as fp:
    json.dump(cam_cal_dict, fp, indent=2)


print("\n\nCamera parameter file saved to: " + str(full_cam_file_path) +"\n\n")






