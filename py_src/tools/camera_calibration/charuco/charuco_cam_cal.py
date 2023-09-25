#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
charuco_cam_cal.py

This program is intended to:
* ingest images of ChArUco targets
* create feature points of the checker corners
* generate camera cal parameters
* save off camera params

''' 

################################
#LOAD LIBRARIES
################################
import os
import cv2
import sys
import glob
import json
import numpy as np


################################
#USER INPUT
################################
show_plots = True
img_extension = '.jpg'

# Define ChArUco target
number_of_target_rows = 5
number_of_target_columns = 5
aruco_square_size = 2.25*0.0254 #m
checker_square_size = 4.5*0.0254 #m


################################
#MAIN CODE
################################
# termination criteria
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

# Defining the dimensions of checkerboard
CHECKERBOARD = (number_of_target_columns-1,number_of_target_rows-1)

# Creating vector to store vectors of 3D points for each checkerboard image
checker_world_points_tot = []
# Creating vector to store vectors of 2D points for each checkerboard image
checker_img_points = [] 
charuco_img_points = []
charuco_ids_tot = []

# Defining the world coordinates for 3D points
checker_world_points = np.zeros((1, CHECKERBOARD[0] * CHECKERBOARD[1], 3), np.float32)

checker_world_points[0,:,:2] = checker_square_size*np.mgrid[0:CHECKERBOARD[0], 0:CHECKERBOARD[1]].T.reshape(-1, 2)
prev_img_shape = None

# Create the ChArUco board dictionary
dictionary = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_1000)
Charuco_board = cv2.aruco.CharucoBoard_create(number_of_target_columns, number_of_target_rows, checker_square_size, aruco_square_size, dictionary)

# Extracting path of individual image stored in a given directory
images = glob.glob('*'+img_extension)

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
found_corner_ctr = 0

for fname in images:

    print("\n" + str(fname))

    img = cv2.imread(fname)

    width = int(img.shape[1] * scale_percent / 100)
    height = int(img.shape[0] * scale_percent / 100)
    dim = (width, height)
    img = cv2.resize(img, dim)

    gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)

    # Detect ArUco markers
    corners, ids, rejectedImgPoints = cv2.aruco.detectMarkers(gray, dictionary)
    # print(corners)

    """
    If desired number of corner are detected, refine the pixel coordinates and display 
    them on the images of checker board
    """
    if ids is not None:
        # Refine the detected corners and obtain the ChArUco corners and ids
        num_corners, charuco_corners, charuco_ids = cv2.aruco.interpolateCornersCharuco(
            corners, ids, gray, board=Charuco_board)

        print("    Detected "+str(len(corners))+" corners")
        # If the entire ChArUco board is detected
        if charuco_corners is not None and num_corners==16:
            found_corner_ctr += 1
            print("    FOUND all ArUco markers in: "+fname)

            charuco_img_points.append(charuco_corners)
            charuco_ids_tot.append(charuco_ids)

            # Extract the checkerboard corners
            checkerboard_corners = []
            for corner in charuco_corners:
                checkerboard_corners.append(corner[0])

            checker_world_points_tot.append(checker_world_points)
            checkerboard_corners = np.asarray(checkerboard_corners)
            corners2 = cv2.cornerSubPix(gray, checkerboard_corners, (11,11),(-1,-1), criteria)
            checker_img_points.append(corners2)
            #print(checkerboard_corners)

            # Draw the detected corners on the image
            if show_plots:
                #img = cv2.aruco.drawDetectedCornersCharuco(img, charuco_corners, charuco_ids)
                img = cv2.drawChessboardCorners(img, CHECKERBOARD, checkerboard_corners, True)

                # Show the image
                cv2.imshow('ChArUco Reprojection', img)
                cv2.waitKey(1)
    else:
        print("    DID NOT FIND corners in: "+fname)

print("\n\n    Done finding corners... attempting calibration! (Using "+str(found_corner_ctr)+" of "+str(len(images))+" total images)")
 
"""
Performing camera calibration by 
passing the value of known 3D points (checker_world_pointsoints)
and corresponding pixel coordinates of the 
detected corners (checker_img_points)
"""
if(charuco_img_points):
    RMS_reproj_err_pix, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(
        checker_world_points_tot, checker_img_points, gray.shape[::-1], None, None,flags=cv2.CALIB_FIX_PRINCIPAL_POINT)

    print("\nReprojection Error (pix, RMS): " + str(RMS_reproj_err_pix))
    print("\nCamera Matrix: "+ str(mtx))
    print("\nDist: " + str(dist))
    print("\nRvecs: " + str(rvecs))
    print("\nTvecs: " + str(tvecs))
    print("")

    image_size = (img.shape[1],img.shape[0]) # tuple required, input args flipped

    #populate dict
    dist_l=dist.tolist()
    newcameramtx_l = mtx.tolist()
    #in this case, cy is vp and cx is up.
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
    charuco_dir = os.path.dirname(os.path.realpath(__file__))
    cam_cal_dir = os.path.dirname(charuco_dir)
    tools_dir = os.path.dirname(cam_cal_dir)
    py_src_dir = os.path.dirname(tools_dir)
    repo_dir = os.path.dirname(py_src_dir)
    cam_config_dir = os.path.join(repo_dir,'data','cam_config')
    full_cam_file_path = os.path.join(cam_config_dir,usr_in+'.json')

    full_cam_file_path=os.path.join(usr_in+'.json')

    with open(full_cam_file_path, 'w') as fp:
        json.dump(cam_cal_dict, fp, indent=2)


print("\n\ncamera parameter file saved to: " + str(full_cam_file_path) +"\n\n")


if show_plots:
    cv2.destroyAllWindows()
