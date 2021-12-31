#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
This code uses a user-defined database, image directory, image file
extension, and darkframe to loop through the Tetra database and generate
inputs for the opencv camera calibration routine
"""

################################
#LOAD LIBRARIES
################################
import os
import sys
import cv2
import json
import time
import pathlib
import numpy as np
from PIL import Image
from scipy.io import savemat
from tetra3_Cal import Tetra3


################################
#USER INPUT
################################
db_name = 'test_db'
path_to_images = ''
image_file_extension = '.jpg'
darkframe_filename = ''
estimated_full_angle_fov = 20

verbose = True
################################
#MAIN CODE
################################

# figure out if a database exists for Tetra
if isinstance(db_name, str):
    path = (pathlib.Path(__file__).parent / db_name).with_suffix('.npz')
else:
    path = pathlib.Path(path).with_suffix('.npz')

# if not, create one
if not os.path.exists(path):
    print("\nUnable to find specified database: " + str(path))
    print("    Generating database...")
    t3 = Tetra3()  #this initializes stuff
    t3.generate_database(max_fov = estimated_full_angle_fov, save_as=db_name)  #this builds a database
    print("    ...complete\n\n")

# init Tetra w/ database
t3 = Tetra3(db_name)

# find darkframe
if darkframe_filename == '':
    print("\nDarkframe file not provided, proceeding without darkframe subtraction.\n")
    darkframe_filename = None
else:
    darkframe_filename = os.path.join(path_to_images, darkframe_filename)
    if os.path.exists(darkframe_filename):
        print('\nDarkframe found!  Applying to calibration images...\n')
    else:
        darkframe_filename = None
        print("\nUnable to find darkframe, proceeding without darkframe subtraction.\n")

# define variables and start processing
path = pathlib.Path(path_to_images)
solution_list = []
AllProjs = np.array([[0,0,0]])
AllCents = np.array([[0,0]])
num_images = 0
for impath in path.glob('*'+image_file_extension):

        num_images+=1
        print('Attempting to solve image: ' + str(impath))
        start_time = time.time()

        img = cv2.imread(str(impath), cv2.IMREAD_GRAYSCALE)
        if img is None:
            print("ERROR ["+str(__name__)+"]:Image file "+impath+" does not exist in path. ")
            sys.exit()
        if verbose: print('Loaded image in {} seconds'.format(time.time()-start_time))
        image_size = (img.shape[1],img.shape[0]) # tuple required, input args flipped

        if darkframe_filename is not None:
            darkframe = cv2.imread(darkframe_filename, cv2.IMREAD_GRAYSCALE)
            img = cv2.subtract(img, darkframe)

        solved = t3.solve_from_image(img, fov_estimate=estimated_full_angle_fov)  # Adding e.g. fov_estimate=11.4, fov_max_error=.1 may improve performance

        try:
            Vecs = []
            ProjVecs = []
            R = solved['Rmat']
            Cents = np.array(solved['MatchCentroids'])
            Cents[:,[0,1]] = Cents[:,[1,0]]
            #print(Cents)  #Swap rows because tetra uses centroids in (y,x)
            AllCents = np.append(AllCents, Cents, axis = 0)
            i=0
            angY = 90 * (np.pi/180)
            RotY = np.array([[np.cos(angY), 0, -np.sin(angY)],
                [0,1,0],
                [np.sin(angY), 0, np.cos(angY)]])
            angZ = -90 * (np.pi/180)
            RotZ = np.array([[np.cos(angZ), -np.sin(angZ), 0 ], 
                [np.sin(angZ), np.cos(angZ), 0],
                [0,0,1]])

            for tup in solved['MatchVectors']:
                v = np.array(tup[1]).transpose()
                #print(v)
                vcam = np.dot(RotZ, np.dot(RotY, np.dot(R, v)))
                #Project Vector onto z = 1 
                proj = np.array([(vcam[0]/vcam[2]),(vcam[1]/vcam[2]), 0])
                #print(proj)
                AllProjs = np.append(AllProjs, [proj], axis = 0)
                ProjVecs.append(proj)
                Vecs.append(vcam)
                i+=1

            imgdict = {'R_inertial_to_camera_guess': R, 'uvd_meas': Cents, 'CAT_FOV':Vecs, 'Projections':ProjVecs, 'NumStars': len(Vecs), 'FOV': solved['FOV']}
            solution_list.append(imgdict)
            print("    ...success! (took " +str(time.time()-start_time)[:8]+" seconds)")

        except:
            print('    ...FAILED to solve \n')

# if some of the images were solved, use the solutions to calibrate the camera
if len(solution_list) > 0:

    print("\n\nFound ("+str(len(solution_list))+") solutions out of (" + str(num_images)+") images\n")

    solved = solution_list[0] # these seem to be fairly consistent, though one day we may want to average all
    fovguess = solved['FOV']
    flpix = image_size[0]/(2*np.tan(np.deg2rad(fovguess)/2))
    matGuess = np.array([[flpix, 0,(image_size[0]/2) - 0.5 ], [0, flpix, (image_size[1]/2)-0.5], [0, 0, 1]])

    AllCents = np.delete(AllCents, 0, axis = 0).astype('float32')
    AllProjs = np.delete(AllProjs, 0, axis = 0).astype('float32')
    RMS_reproj_err_pix, mtx, dist, rvecs, tvecs = cv2.calibrateCamera([AllProjs], [AllCents], image_size, matGuess, None, flags=(cv2.CALIB_USE_INTRINSIC_GUESS+cv2.CALIB_FIX_PRINCIPAL_POINT))
    [fovx, fovy, fl, pp, ar] = cv2.calibrationMatrixValues(mtx, image_size, 4.8e-3*image_size[0], 4.8e-3*image_size[1] ) # Note that this may not work without an accurate detector size (Args 3 and 4)

    print("\nReprojection Error (pix, RMS): " + str(RMS_reproj_err_pix))
    print("\nFoV (x,y): "+str(fovx)+", "+str(fovy))
    print("Focal Length (pix): " + str(flpix))
    print("Focal Length (mm): " + str(fl))
    print("\nCamera Matrix: "+ str(mtx))
    print("\nDist: " + str(dist))
    print("\nRvecs: " + str(rvecs))
    print("\nTvecs: " + str(tvecs))
    print("")

    # populate dict
    dist_l=dist.tolist()
    newcameramtx_l = mtx.tolist()
    # in this case, cy is vp and cx is up.
    cam_cal_dict = {'camera_matrix': newcameramtx_l, 'dist_coefs': dist_l, 'resolution':[image_size[0],image_size[1]], 'camera_model':'Brown','k1':dist_l[0][0],'k2':dist_l[0][1],'k3':dist_l[0][4],'p1':dist_l[0][2],'p2':dist_l[0][3],'fx':newcameramtx_l[0][0],'fy':newcameramtx_l[1][1],'up':newcameramtx_l[0][2],'vp':newcameramtx_l[1][2],'skew':newcameramtx_l[0][1], 'RMS_reproj_err_pix':RMS_reproj_err_pix}

    # save
    usr_in = 'generic_cam_params'
    usr_in_split = usr_in.split('.json')
    usr_in = usr_in_split[0]
    cam_config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))),'data','cam_config')
    full_cam_file_path = os.path.join(cam_config_dir,usr_in+'.json')
    with open(full_cam_file_path, 'w') as fp:
        json.dump(cam_cal_dict, fp, indent=2)

    print("\n\nCamera parameter file saved to: " + str(full_cam_file_path) +"\n\n")



else: print("\n\n\nNo solutions found (at all).  Exiting unsuccessfully...\n\n\n")
