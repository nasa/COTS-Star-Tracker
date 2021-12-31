#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
support_functions.py

File containing functions used for a variety of tasks
'''

################################
#SUPPORT FUNCTIONS
################################
def timing_decorator(func):
    import functools
    import time
    @functools.wraps(func)
    def timing_wrapper(*args, **kwargs):
        t_initial = time.time()
        output = func(*args, **kwargs)
        t_end = time.time()
        print("    [TIME DECORATOR]: the ({}) function completed in {} seconds".format(func.__name__, round((t_end - t_initial),5)))
        return output
    return timing_wrapper


def stars_in_fov(u_i, t_i2c, c, nrow, ncol, fov):
    import math
    import numpy as np
    import array_transformations as xforms
    u_c = xforms.vector_array_transform(t_i2c, u_i)
    # Get stars within circular FOV
    idx_uc_fov = np.nonzero(u_c[:][2] >= math.cos(fov))[0]
    u_fov_c = u_c[:, idx_uc_fov]

    # Convert homogeneous coordinates to pixel coordinates
    p = xforms.vector2pixel(c, u_fov_c)

    # Get indices of pixels in range of rectangular fov
    idx_rect_fov = np.where(np.logical_and(
        np.logical_and(
            p[0, :] <= ncol, p[0, :] >= 1),
        np.logical_and(
            p[1, :] <= nrow, p[1, :] >= 1)))[0]

    # Grab original indices from list
    star_idx_fov = idx_uc_fov[idx_rect_fov]

    return star_idx_fov


def basic_cam(): #this function returns some arbitrary camera parameters as a stand-in when real values aren't available to support troubleshooting
    import cam_matrix as cm
    # Test camera option 1
    # ncol = 640
    # nrow = 480
    # dx = 583.2829786373293
    # up = 320.0  # up = (ncol+1)/2;
    # dy = 579.4112549695428
    # vp = 240.0  # vp = (nrow+1)/2;
    # sk = 0.0
    # fx = 7.891092810360696e+03
    # fy = 7.906007692160952e+03

    # Test camera option 2
    nrow = 2048
    ncol = 2592
    dx = 7.891092810360696e+03
    dy = 7.906007692160952e+03
    up = 1.274779317743486e+03
    vp = 1.066692381188306e+03
    sk = 0
    c = cm.cam_matrix(dx, dy, up, vp, sk)
    cinv = cm.cam_matrix_inv(c)
    fov = cm.cam2fov(cinv, nrow, ncol)
    return c, nrow, ncol, fov


def catalog2fov(q, c, img, stars):
    import cam_matrix as cam
    import array_transformations as xforms
    import numpy as np
    cam_res = img.shape
    nrow = cam_res[0]
    ncol = cam_res[1]
    cinv = cam.cam_matrix_inv(c)
    fov = cam.cam2fov(cinv, cam_res[0], cam_res[1])
    t_hat = xforms.quat2attitude_matrix(q[0:3], q[3])
    x_cat_t = xforms.vector_array_transform(t_hat, stars)

    x_cat_fov_idx = np.where(x_cat_t[2, :] >= np.cos(fov))[0]
    x_cat_fov = x_cat_t[:, x_cat_fov_idx]
    p_fov = xforms.vector2homogeneous(c, x_cat_fov)
    idx_rect = np.where((p_fov[0, :] <= ncol) & (p_fov[1, :] <= nrow) &
                        (p_fov[0, :] > 0) & (p_fov[1, :] > 0))[0]
    return p_fov[:, idx_rect]


def reproject(img, c, idmatch, q_est, x_obs=None, x_cat=None):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import array_transformations as xforms
    import cam_matrix as cam
    import numpy as np

    c_inv = cam.cam_matrix_inv(c)
    fov = cam.cam2fov(c_inv,img.shape[0],img.shape[1])

    id_cat = np.where(idmatch != -1)[0]
    nid_cat = np.where(idmatch == -1)[0]
    x_cat_found = x_cat[:, idmatch[id_cat][:, 0]]

    mask = np.ones(len(x_cat[0]), dtype=bool)
    mask[idmatch[id_cat][:, 0]] = False
    x_cat_nfound = x_cat[..., mask]

    cat_not_found = catalog2fov(q_est, c, img, x_cat_nfound)
    cat_found = catalog2fov(q_est, c, img, x_cat_found)

    x_obs_found = x_obs[:, id_cat]
    x_obs_invalid = x_obs[:, nid_cat]

    dpi = 100
    height, width = img.shape
    figsize = width / float(dpi), height / float(dpi)

    # Create a figure of the right size with one axes that takes up the full figure
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.autoscale(tight=True)
    # plt.axis('off')
    imgplot = plt.imshow(img, cmap='Greys_r')

    row = 1
    col = 0
    if x_obs_found is not None:
        p_obs = xforms.vector2homogeneous(c, x_obs_found)
        plt.plot(p_obs[col, :], p_obs[row, :], 'cs',
                 markersize=30, markerfacecolor="None",
                 label="Matched candidates")
    if x_obs_invalid is not None:
        p_inv = xforms.vector2homogeneous(c, x_obs_invalid)
        plt.plot(p_inv[col, :], p_inv[row, :], 'r<',
                 markersize=20, markerfacecolor="None",
                 label="Invalid candidates")
    if cat_found is not None:
        plt.plot(cat_found[col, :], cat_found[row, :], 'gd',
                 markersize=25, markerfacecolor="None",
                 label="Matched catalog stars")
    if cat_not_found is not None:
        plt.plot(cat_not_found[col, :], cat_not_found[row, :], 'ro',
                 markersize=25, markerfacecolor="None",
                 label="Unmatched catalog stars")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='lower left', prop={'size': 20})
    plt.show()


def find_candidate_stars(img, low_thresh, hi_thresh,
                         min_star_area, max_star_area, graphics=False):
    import numpy as np
    import cv2
    from rpi_core import calculate_center_intensity

    if low_thresh is None:
        percent_thresh = 0.995
        hist, bins = np.histogram(img.flatten(), 256, [0, 256])
        cdf = hist.cumsum()
        cdf_normalized = cdf / cdf.max()
        idx_first_bin = np.where(cdf_normalized>percent_thresh)[0][0]
        low_thresh = idx_first_bin

    if hi_thresh is None:
        dtype = img.dtype
        hi_thresh = np.iinfo(dtype).max

    # convert the grayscale image to binary image
    _, im_binary = cv2.threshold(img, low_thresh, hi_thresh, cv2.THRESH_BINARY)

    # find contours in the binary image
    opencv_ver = cv2.__version__
    if int(opencv_ver[0]) == 3:
        junk_img, contours, hierarchy = cv2.findContours(
           im_binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE) #opencv3
    else:
        contours, hierarchy = cv2.findContours(
            im_binary, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE) #opencv2 and 4

    if graphics:
        import matplotlib.pyplot as plt
        #cv2.drawContours(img, contours, -1, (0,255,0), 3) TODO look into fixing this
        #cv2.imshow('wow1',img) #this did not result in contours drawing.  Would be nice to somehow show contours
        plt.imshow(im_binary), plt.show()


    print("Found " + str(len(contours)) + " contours")

    try:
        connectivity = 8
        _, _, stats, _ = cv2.connectedComponentsWithStats(
            im_binary, connectivity, cv2.CV_32S)
        # Set the background region area to zero area
        stats[0, cv2.CC_STAT_AREA] = 0

        # TODO determine if this needs to be resurrected
        # TODO should this use the raw image or binary image?  I think it should use raw like it is...
        coi, intensities = calculate_center_intensity(img, stats, min_star_area, max_star_area)
        print("Found " + str(len(coi)) + " centroids")
        return coi, intensities


    except Exception as exc:
        raise Exception("ERROR ["+str(__name__)+"]: An unknown error occurred in {0}. {1}".format(__name__, exc.args[-1]))

    '''
    finally: #TODO: should we remove this "finally" gate and keep the indented contents?
        # centroid blobs
        centroids = []
        intensity = []
        for c in contours:
            # calculate moments for each contour
            M = cv2.moments(c)

            # calculate x,y pixel location of center
            if M["m00"] > 0:
                cx = M["m10"] / M["m00"]
                cy = M["m01"] / M["m00"]
                centroids += [[cx, cy]]
                intensity += [M["m00"]]

        centroids = np.array(centroids)
        intensity = np.array(intensity)
        print("Found " + str(len(centroids)) + " centroids")
        #print(centroids)

        return centroids, intensity
     '''




def deg2dms(angle, in_str=None, img_filename=None, conversion=1):
    arcsec = angle*conversion
    if in_str is None:
        in_str = 'Quaternion Error'
    if img_filename is not None:
        in_str = img_filename + ': ' + in_str
    # print(in_str, str(deg), "degrees,", str(arcmin),
    #           "',", str(arcsec), "''")
    print(in_str, str(arcsec), " arcmin")
    return


def att2angle(R):
    import numpy as np
    import scipy.linalg as linalg
    theta = np.arccos((np.trace(R)-1)/2)
    t_arr = np.array([
        [R[1, 2] - R[2, 1]],
        [R[2, 0] - R[0, 2]],
        [R[0, 1] - R[1, 0]]])
    e = 1/(2*np.sin(theta))*t_arr
    phi = theta*e
    phix = np.array([
        [0, float(-phi[2]), float(phi[1])],
        [float(phi[2]), 0, float(-phi[0])],
        [float(-phi[1]), float(phi[0]), 0]])
    t = linalg.expm(-phix)

    # l = linalg.logm(R.transpose(), True)
    # phi = np.array([[l[2, 1]],[l[0, 2]],[l[1, 0]]])
    # phix = np.array([[0, -phi[2], phi[1]],[phi[2], 0, -phi[0]],[-phi[1], phi[0], 0]])
    # t = linalg.expm(-phix)
    return phi, theta

def attitude_error(R1, R2):
    # Using formula from http://www.boris-belousov.net/2016/12/01/quat-dist/
    import numpy as np
    R = np.dot(R1, R2.T)
    # x, y, z, theta = att2angle(R)
    return att2angle(R)


def quat_error(p, q):
    # Using formula from http://www.boris-belousov.net/2016/12/01/quat-dist/
    import numpy as np
    R = np.dot(p, q.T)
    return arccos((np.trace(R)-1)/2)
