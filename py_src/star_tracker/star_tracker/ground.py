#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
ground.py

File containing functions used on ground

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

'''

################################
#LOAD LIBRARIES
################################
from array_transformations import check_axis_decorator

################################
#SUPPORT FUNCTIONS
################################
def kvector(input_cat):
    import numpy as np
    import rpi_core
    m, q, sorted_cat = rpi_core.kvec_values(input_cat)
    nrow = len(sorted_cat)

    k = np.ones((nrow, 1), dtype=int)

    k[0] = 0
    k[-1] = nrow
    for x in np.arange(1, nrow-1, 1, dtype=int):
        # Create k-vector catalog by finding the number of
        # items in the original catalog that are below the
        # value of the line
        l = k[x-1][0]  # grab the previous (smaller) element in the array
        z = m*(x+1) + q  # eqn 1
        for y in np.arange(k[x-1], nrow-1, 1, dtype=int):

            # If the calculated z matches/exceeds that of the current catalog
            # entry, increment l by 1
            if z >= sorted_cat[y, 1]:  # (eqn 2)
                l += 1
            else:
                break
        k[x] = l
    return k, m, q, sorted_cat


@check_axis_decorator(2)
def equatorial2vector(ra_de, axis=None):
    # Axis indicates if array is 2xn array (axis=0), or nx2 array (axis=1)
    # Axis = 0: each column represents individual equatorial elements (RA or DE)
    # Axis = 1: each row represents individual equatorial elements (RA or DE)
    import numpy as np

    ra_de_rad = np.radians(ra_de)
    if ra_de_rad.ndim == 1:
        ra = np.array(ra_de_rad[0])
        de = np.array(ra_de_rad[1])
    else:
        if axis == 0:
            # axis = 0 2x3
            #  | RA RA RA |
            #  | DE DE DE |
            ra = ra_de_rad[0, :]
            de = ra_de_rad[1, :]
        elif axis == 1:
            # axis = 1 4x2
            #  | RA DE |
            #  | RA DE |
            #  | RA DE |
            #  | RA DE |
            ra = ra_de_rad[:, 0]
            de = ra_de_rad[:, 1]
        else:
            raise
    los_vector = np.array([np.cos(de)*np.cos(ra),
                           np.cos(de)*np.sin(ra),
                           np.sin(de)])
    return los_vector, ra, de


@check_axis_decorator(3)
def vector2equatorial(v, axis=None):
    # Axis indicates if array is 3xn array (axis=0), or nx3 array (axis=1)
    # Axis = 0: each star's data is a column
    # Axis = 1: each star's data is a row
    import numpy as np

    if axis == 0:
        de = np.arcsin(v[2, :])
        ra = np.arctan2(v[1, :], v[0, :])
        ra_de = np.vstack((ra, de))
    elif axis == 1:
        de = np.arcsin(v[:, 2])
        ra = np.arctan2(v[:, 1], v[:, 0])
        ra_de = np.hstack((ra, de))
    return np.degrees(ra_de)


@check_axis_decorator(2)
def lpq_orthonormal_basis(ra_de_deg, axis=None):
    # Axis indicates if array is nx2 array (axis=0), or 2xn array (axis=1)
    # Axis = 0: each star's data is a column
    # Axis = 1: each star's data is a row
    import numpy as np

    los, ra_rad, de_rad = equatorial2vector(ra_de_deg, axis=axis)

    p_hat = np.array([
        -np.sin(ra_rad),
        np.cos(ra_rad),
        np.zeros(len(ra_rad))])
    q_hat = np.array([
        -np.sin(de_rad) * np.cos(ra_rad),
        -np.sin(de_rad) * np.sin(ra_rad),
        np.cos(de_rad)])
    return los, p_hat, q_hat


@check_axis_decorator(2)
def proper_motion_correction(ra_de_deg, pm_rade_mas, plx_mas, rB, t, t_ep, axis=None):
    import array_transformations as xforms
    import math
    import numpy as np
    import numpy.matlib as matlib

    if ra_de_deg.shape != pm_rade_mas.shape:
        raise ValueError("ERROR ["+str(__name__)+"]: Dimensions for RA/DEC and proper motion do not agree")

    deg2rad = math.pi/180;
    arcsec2deg = 1/3600;
    mas2arcsec = 1/1000;
    mas2rad = mas2arcsec*arcsec2deg*deg2rad;

    AU = 1.496e8; #in km
    # return 3xn Line-of-Sight vector array and orthogonal
    # unit vectors p and q so that {p,q,l} form a
    # right-handed orthonormal basis
    los, p_hat, q_hat = lpq_orthonormal_basis(ra_de_deg, axis=axis)

    # Convert proper motion of right ascension and declination to rad/yr 2xn array
    pm_rade_rad = mas2rad*pm_rade_mas
    mu_ra = pm_rade_rad[0, :]
    mu_de = pm_rade_rad[1, :]
    pm = (t - t_ep)*(mu_ra*p_hat + mu_de*q_hat)

    # Convert parallax to rads 1xn array
    plx_rad = mas2rad*plx_mas

    # Define au (Astronomical Units) in kilometers
    au = 1.496e8  # km

    # Find location of observer in units of AU
    if rB is None:
        rB = np.array([[149597870.693],[0],[0]])
    rObs_AU = rB/au;

    # Define rB as BCRF position (in km) of celestial object that the spacecraft orbits
    # if no rB provided, assume that s/c is orbiting Earth: 149597870.693 km

    # Incorporate proper motion into already existing u vector:
    # ui = l_i+(t-t_ep)*(mua_i*p+mud_i*q)-(w_i*rB)/AU
    r_au_mat = matlib.repmat(rObs_AU, 1, len(plx_rad))
    plx = -plx_rad*r_au_mat;

    u = los + pm + plx
    return xforms.normalize_vector_array(u), los


def get_Hipparcos_data(starcat_file, magnitude_name, excess_rows,
                       brightness_thresh=None, index_col=2):
    import pandas as pd
    # read catalog file using Pandas
    starcat = pd.read_csv(starcat_file,skiprows=excess_rows,sep='\t',
                          comment='#',index_col=index_col)
    # return all stars brighter than given brightness threshold
    if brightness_thresh is None:
        return starcat
    else:
        return starcat[starcat[magnitude_name]<brightness_thresh]


def read_star_catalog(starcat_file, brightness_thresh, t=None, cat_ep=None,
                      rB=None, excess_rows=None, index_col=2):
    from astropy.time import Time
    import numpy as np
    import ground

    if t is None:
        t = Time.now().byear
    else:
        t = Time(t).byear
    # t = 2024.25  TODO: why is this here?  What was this used for?
    if cat_ep is None:
        cat_ep = 1991.250
    cat_ep = Time(cat_ep, format='byear').byear

    # Column header names in Hipparcos star catalog
    HIP_ID = 'HIP'
    MAG = 'Hpmag'
    # RAJ2000 = '_RAJ2000'
    # DEJ2000 = '_DEJ2000'
    RA = 'RArad'
    DE = 'DErad'
    PLX = 'Plx'
    PM_RA = 'pmRA'
    PM_DE = 'pmDE'
    # B_V = 'B-V'

    starcat = get_Hipparcos_data(starcat_file, MAG, excess_rows,
                                 brightness_thresh=brightness_thresh,
                                 index_col=index_col)

    ra_de = np.array([starcat[RA].values, starcat[DE].values])
    pm_rade = np.array([starcat[PM_RA].values, starcat[PM_DE].values])
    plx = np.array(starcat[PLX].values)

    # Correct catalog entries for proper motion: 'u' is array of unit vectors to each star
    u, _ = ground.proper_motion_correction(ra_de, pm_rade, plx, rB, t, cat_ep)
    return u, starcat


def create_star_catalog(starcat_file, brightness_thresh, cat_ep=None, t=None, rB=None,
                        excess_rows=None, index_col=2, save_vals=False, fov=None, save_dir=None):
    import os
    import rpi_core
    import numpy as np
    import itertools as it

    # Correct catalog entries for proper motion
    u, _ = read_star_catalog(
        starcat_file, brightness_thresh, excess_rows=excess_rows,
        cat_ep=cat_ep, t=t, rB=rB, index_col=index_col)

    # Create star pairs using nchoosek
    # star_idx = np.arange(0, len(u[0]))
    star_idx = np.arange(0, u.shape[1])
    star_pairs = np.array( list(it.combinations(star_idx, 2)) )

    # Form star pairs unit vectors into an nx6 array
    u_starpairs = np.vstack((u[:, star_pairs[:, 0]],
                             u[:, star_pairs[:, 1]]))

    # Calculate interstar angles
    istar_angle = rpi_core.interstar_angle(u_starpairs)

    # Remove star pairs that fall outside field of view angle
    if fov is not None:
        sp_fov = np.where(istar_angle < 2.0*fov)[0]
        istar_angle = istar_angle[sp_fov]
        star_pairs = star_pairs[sp_fov, :]
        del sp_fov, fov

    # Create star pair catalog with interstar angle
    k, m, q, isa_cat = kvector(istar_angle)

    isa_cat_idx = np.hstack((star_pairs, isa_cat))

    if save_vals:
        if save_dir == None:
            print("[CREATE_STAR_CATALOG]: no directory provided, exiting without saving")
            return k, m, q, u, isa_cat_idx

        # Save values or return them if in a function
        np.save(os.path.join(save_dir,'k'), k)
        np.save(os.path.join(save_dir,'m'), m)
        np.save(os.path.join(save_dir,'q'), q)
        np.save(os.path.join(save_dir,'u'), u)
        # np.save(os.path.join(save_dir,'isa_cat'), isa_cat)
        np.save(os.path.join(save_dir,'indexed_star_pairs'), isa_cat_idx)

    return k, m, q, u, isa_cat_idx


def create_darkframe(img_list, numImages):
    import cv2 as cv
    import numpy as np
    # TODO: use GLOB to read images better
    # TODO: explicitly use filenames instead of grabbing the first n images?
    # TODO: use context manager to open images

    n_images = len(img_list)
    if n_images < numImages: print("Length of image list {0} is less than requested "
          "number of images {1}".format(n_images, numImages))

    numImages = min(numImages, n_images)
    if len(img_list) > numImages: img_list = img_list[0:numImages]
    if numImages == 0:
        print('No images found for dark frame creation.')
        return None

    img_sample = cv.imread(img_list[0], cv.IMREAD_GRAYSCALE)
    img_size = img_sample.shape[0:2]
    img_dtype = img_sample.dtype
    images = np.zeros((img_size[0], img_size[1], numImages))

    for count, filename in enumerate(img_list):
        img = cv.imread(filename, cv.IMREAD_GRAYSCALE)
        if img == None:
            raise('Image {0} (count={1}) is None'.format(filename, count))
        images[:, :, count] = img
            # count += 1

    if images == None:
        print('No images found for darkframe generator. Returning None')
        return None
    images_median = np.median(images, axis=2)
    return images_median.astype(img_dtype)


def create_distortion_map(camera_json, distortion_map_path, save_dist_map=True):
    import numpy as np
    import cam_matrix
    import array_transformations as xforms
    # this function only allows Brown distortion at the moment
    c, img_size, dist_coefs = cam_matrix.read_cam_json(camera_json)
    c_inv = cam_matrix.cam_matrix_inv(c)
    # img_size must be the reverse of that listed in MATLAB for sub2ind and ind2sub to play nice
    # img_size = array([img_size[1], img_size[0]])

    #dist_coefs = array([k_1, k_2, p_1, p_2, k_3])
    k1 = dist_coefs[0]
    k2 = dist_coefs[1]
    k3 = dist_coefs[4]
    p1 = dist_coefs[2]
    p2 = dist_coefs[3]

    # cam_resolution, [ncols, nrows]
    nrows = img_size[1]
    ncols = img_size[0]
    # img_size = flip(img_size)

    # Initalize holding arrays for interpolation variables
    idxUndistorted = np.zeros([nrows*ncols, 1], dtype=int)
    idxDistorted   = np.zeros([nrows*ncols, 4], dtype=int)
    interpWeights  = np.zeros([nrows*ncols, 4])

    # Loop through every pixel in undistorted image and obtain bilinear interpolation weights
    # Bilinear interpolation follows Eqs. 3.6.1-3.6.5 from [Press et al, 2007]

    #Set tolerance for repeated pixels
    #subpix_tol = 0.0001;

    #Initalize pixel counter
    # DO NOT use 'enumerate' in this case since count requires that certain
    # parameters are met to increment
    count = 0;

    # Using 0-indexing
    indexing = 0
    for row in np.arange(indexing, nrows+indexing):
        for col in np.arange(indexing, ncols+indexing):
            # Get current pixel coordinates in undistorted image.
            # Written in homogeneous coordinates.
            uv = np.array( [[col], [row]] ) + 1-indexing
            uv_h = xforms.camera2homogeneous(uv)
            # uv1 = array([[col], [row], [1.0]])

            # Convert undistorted pixel coord to xy coord
            # The inverse of Eq. 7 in [Christian et al, 2016].
            xy = xforms.vector_array_transform(c_inv, uv_h)

            # Isolate x and y coordinates in undistorted image
            x = xy[0, 0];
            y = xy[1, 0];

            # Compute r^2 in undistorted coordinates
            r2 = x**2 + y**2;

            # Additional even power of r needed for Brown model
            r4 = r2*r2;
            r6 = r4*r2;

            # See pp 375-376 of [Bradski & Kaehler, 2008].
            # Also see Eq. 6 of [Christian et al, 2016].
            xy_dist_pre = (1 + k1*r2 + k2*r4 + k3*r6) * np.array([x, y])
            xy_distorted_post = np.array([2*p1*x*y + p2*(r2+2*x**2), p1*(r2+2*y**2) + 2*p2*x*y])
            xy_distorted = xforms.camera2homogeneous(xy_dist_pre + xy_distorted_post)

            # Get distorted uv coordinate
            # See Eq. 7 of [Christian et al, 2016].
            uvDistorted = xforms.vector_array_transform(c, xy_distorted)-1


            # Find indices for the four pixels that surround the query pixel the raw (distorted) image
            # Remember that [1,1] is upper lefthand corner of image.
            # Column number is u-direction. Row number is v-direction.
            colLeft  = np.floor(uvDistorted[0]).astype('int') # From [Press et al, 2007]: x_{1i}
            colRight = np.ceil( uvDistorted[0]).astype('int') # From [Press et al, 2007]: x_{1(i+1)}
            rowUp    = np.floor(uvDistorted[1]).astype('int') # From [Press et al, 2007]: x_{2j}
            rowDown  = np.ceil( uvDistorted[1]).astype('int') # From [Press et al, 2007]: x_{2(j+1)}

            # Check to make sure all four of the query pixels in the new image
            # lie inside of the original (distorted) image.
            if rowUp >= 0 and colLeft >= 0 and rowDown < nrows and colRight < ncols:

                # Store index of current pixel
                order = 'F'
                output_fmt = 'P'

                idxUndistorted[count, 0] = xforms.sub2ind(
                    img_size, row, col, order=order, output_fmt=output_fmt)

                # Get the indices of corresponding surrounding pixels
                idxUpperLeft  = xforms.sub2ind(
                    img_size,rowUp,colLeft,order=order, indexing=indexing, output_fmt=output_fmt)
                idxUpperRight = xforms.sub2ind(
                    img_size,rowUp,colRight, order=order, output_fmt=output_fmt)
                idxLowerRight = xforms.sub2ind(
                    img_size,rowDown,colRight, order=order, output_fmt=output_fmt)
                idxLowerLeft  = xforms.sub2ind(
                    img_size,rowDown,colLeft, order=order, output_fmt=output_fmt)

                # Store indicies from distorted image in clockwise order, starting from upper left
                idxDistorted[count, :] = [idxUpperLeft,
                                          idxUpperRight,
                                          idxLowerRight,
                                          idxLowerLeft]

                # Compute distance from distorded point to surrounding pixel centers
                # See Eq. 3.6.4 from [Press et al, 2007]. Note there is no
                # denominator since the distance between adjacent pixels is
                # always unity (we devide by 1). Baseline equation is augmented
                # with check to protect against case where distorted point lies
                # exactly on an integer pixel value. This leads to floor() and
                # ceil() producing the same value and must be handled differently.
                if rowDown==rowUp: # abs(colRight-colLeft)<subpix_tol
                    s = 0.5
                else:
                    # Note: What we call "s" here is called "u" in [Press et al,
                    # 2007]. We choose to use "s" to avoid possibility of
                    # confusion with pixel [u,v] coordinates.
                    s  = uvDistorted[1] - rowUp
                if colRight==colLeft:  # abs(colRight-colLeft)<subpix_tol
                    t = 0.5
                else:
                    t  = uvDistorted[0] - colLeft

                # Get weights based on distances. What we call the interpolation
                # "weights" are coefficients in Eq. 3.6.5 from [Press et al, 2007].
                wtUpperLeft  = (1-t)*(1-s)
                wtUpperRight = t*(1-s)
                wtLowerRight = t*s
                wtLowerLeft  = (1-t)*s

                # Store indicies from distorted image in clockwise order, starting from upper left
                # This ordering scheme is consistent with index storage for idxDistorted
                interpWeights[count,:] = [ wtUpperLeft, wtUpperRight, wtLowerRight, wtLowerLeft ]

                # Increment counter
                count += 1

    # Keep only the usable elements
    idxUndistorted = idxUndistorted[0:count, :]
    idxDistorted = idxDistorted[0:count, :]
    interpWeights = interpWeights[0:count, :]

    np.savez(distortion_map_path, c=c, cam_resolution=img_size, dist_coefs=dist_coefs, idxUndistorted=idxUndistorted,idxDistorted=idxDistorted,interpWeights=interpWeights)
    return c, img_size, dist_coefs, idxUndistorted, idxDistorted, interpWeights


def find_file(name, path):
    # finds file {name} in the {path} using the os.walk command
    import os
    # find file name in path
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)
    else:
        raise FileNotFoundError("ERROR ["+str(__name__)+"]: File {0} not found.".format(name))

def find_files_pattern(pattern, path, exclude=None):
    import os, fnmatch
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))

    if exclude is not None:
        result = [i for i in result if exclude not in i]
    return result
