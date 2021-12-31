#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
main.py

This script is intended to contain the primary function(s)
required to get an attitude estimate from an image
of stars
'''

################################
#LOAD LIBRARIES
################################
from support_functions import timing_decorator, reproject

################################
#SUPPORT FUNCTIONS
################################
@timing_decorator
def star_tracker(img_file_name, cam_config_file_name, m=None, q=None, x_cat=None, k=None,
                 indexed_star_pairs=None, darkframe_file=None, undist_img_bool=True, n_stars=30,
                 low_thresh_pxl_intensity=None, hi_thresh_pxl_intensity=None, min_star_area=4,
                 max_star_area=36, isa_thresh=0.0008, nmatch=6, watchdog=None, graphics=False, verbose=False):

    import os
    import cv2
    import sys
    import time
    import rpi_core
    import numpy as np
    import cam_matrix as cam
    import support_functions
    import array_transformations as xforms

    start_time = time.time()
    # opencv will return None object if file not found
    img = cv2.imread(img_file_name, cv2.IMREAD_GRAYSCALE)

    if img is None:
        print("ERROR ["+str(__name__)+"]:Image file "+img_file_name+" does not exist in path. ")
        sys.exit()
    if verbose: print('\nLoaded image in {} seconds'.format(time.time()-start_time))

    camera_matrix, cam_resolution, dist_coefs = cam.read_cam_json(cam_config_file_name)
    dist_coefs = np.array(dist_coefs)
    im_resolution = np.array([img.shape[1], img.shape[0]], dtype=int) #pixels

    if np.any(im_resolution != cam_resolution):
        raise ValueError('ERROR ['+str(__name__)+']: Input image and camera file do not contain the same camera parameters.  Image: '+str(im_resolution)+', Camera file: '+str(cam_resolution))

    if darkframe_file is not None:
        darkframe = cv2.imread(darkframe_file, cv2.IMREAD_GRAYSCALE)
        img = cv2.subtract(img, darkframe)

    # # undistort points #TODO determine location of this vs CV2 and the functionality difference
    # if undist_img_bool is True:
    #     img = rpi_core.undistort_image(
    #         img, idxUndistorted, idxDistorted, interpWeights)

    if min_star_area >= max_star_area:
        min_star_area = 4
        max_star_area = 36
        if verbose: print('WARNING ['+str(__name__)+']: min_star_area greater than or equal to max_star_area.  Both have been reset to their default values.')


    centroids, intensities = support_functions.find_candidate_stars(
        img, min_star_area=min_star_area, max_star_area=max_star_area,
        low_thresh=low_thresh_pxl_intensity, hi_thresh=hi_thresh_pxl_intensity, graphics=graphics)

    # if fewer than 3 centroids are found, don't bother.  Gums things up downstream.
    if len(centroids) < 3:
        x_obs = np.array(np.nan)
        q_est = np.array(np.nan)
        idmatch = [0]
        nmatches = 0
        if verbose: print("Found too few centroids (< 3) to continue.  Exiting...")
        return q_est, idmatch, nmatches, x_obs, img

    if undist_img_bool is True:
        centroidsUnd = centroids.reshape(len(centroids), 1, 2)
        undistorted_centroids_offsets = cv2.undistortPoints(
            centroidsUnd, camera_matrix, dist_coefs)
        centroids = (centroidsUnd + undistorted_centroids_offsets)
        centroids = centroids.reshape((len(centroidsUnd), 2))

    if graphics is True:
        import matplotlib.pyplot as plt
        mask = np.zeros((img.shape[0],img.shape[1]))
        for centroid in centroids:
            cv2.circle(mask, (int(centroid[0]), int(centroid[1])), 30, (255, 255, 255), 5)
            cv2.circle(mask, (int(centroid[0]), int(centroid[1])), radius=0, color=(255, 255, 255), thickness=-1)

        mask = np.ma.masked_where(mask < 1, mask)
        plt.imshow(img, cmap="Greys_r")
        plt.imshow(mask, cmap="gist_rainbow")
        plt.show()

    intensities = intensities.squeeze()
    original_num_centroids = len(intensities)
    n_stars = min(len(intensities), n_stars)
    bright_star_idx = intensities.argsort()[-n_stars:]
    centroids = centroids[bright_star_idx, :]
    if verbose: print("Using " + str(len(centroids)) + " centroids")


    if graphics is True and (original_num_centroids > n_stars):
        import matplotlib.pyplot as plt
        mask = np.zeros((img.shape[0],img.shape[1]))
        for centroid in centroids:
            cv2.circle(mask, (int(centroid[0]), int(centroid[1])), 30, (255, 255, 255), 5)
            cv2.circle(mask, (int(centroid[0]), int(centroid[1])), radius=0, color=(255, 255, 255), thickness=-1)

        mask = np.ma.masked_where(mask < 1, mask)
        plt.imshow(img, cmap="Greys_r")
        plt.imshow(mask, cmap="gist_rainbow")
        plt.show()

    cinv = cam.cam_matrix_inv(camera_matrix)

    # Corrected centroids in unit vectors
    x_obs = 0
    x_obs = xforms.pixel2vector(cinv, centroids.T)
    if verbose:
        print("x_obs[0] is length "+str(len(x_obs[0])))

    k_vector_interp = (m, q)
    test = np.sum(x_obs)

    if np.isnan(test): 
        q_est = np.array(np.nan)
        idmatch = [0]
        nmatches = 0
        if verbose: print("Failure to calculate x_obs, exiting...")
    else:
        if verbose: print("Starting star ID")

        q_est, idmatch, nmatches = rpi_core.triangle_isa_id(x_obs, x_cat, indexed_star_pairs, isa_thresh,
                                                            nmatch, k, k_vector_interp, watchdog=watchdog, verbose=verbose)

        if graphics is True and nmatches > 0:

            #reproject(img, camera_matrix, x_obs=x_obs_found, x_inv=x_obs_invalid, cat_found=cat_found, cat_nfound=cat_nfound)
            reproject(img, camera_matrix, idmatch, q_est, x_obs=x_obs, x_cat=x_cat)
            # matched candidates: pixel position of the blob that we think is a star
            # invalid candidates: pixel position of a blob we found, but that we don't think are stars
            # matched catalog stars: this is where the star should be assuming the quat is true
            # unmatched catalog stars: this is where a catalog star should be, but we couldn't find it (or maybe it wasn't bright enough to be considered)



    # Return x_obs and the img
    return q_est, idmatch, nmatches, x_obs, img


