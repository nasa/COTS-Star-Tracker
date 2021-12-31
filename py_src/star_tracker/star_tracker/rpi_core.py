#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
rpi_core.py

File containing functions used in flight

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
from support_functions import timing_decorator

################################
#USER INPUT
################################
test_bool = True

################################
#SUPPORT FUNCTIONS
################################
@timing_decorator
def nchoosek(n, k):
    if n < 0 or k < 0:
        return 0
    # Calculates the number of combinations by choosing k number of
    # combinations out of n possible numbers
    # Inspired by Nas Banov's answer in: https://stackoverflow.com/questions/3025162/statistics-combinations-in-python
    import operator as op    # or mul=lambda x,y:x*y
    import fractions, functools
    return int(functools.reduce(op.mul, (fractions.Fraction(n-i, i+1) for i in range(k)), 1))


@check_axis_decorator(6)
def interstar_angle(star_pair, axis=None):
    # Calculates the interstar angle given unit vectors pointing to respective
    # stars
    # Accepts both 3xn (preferred) and nx3 arrays of vectors
    # axis = 0: pairs of stars are column vectors
    # axis = 1: pairs of stars are row vectors
    import numpy as np
    import array_transformations as xforms
    # determine if star pairs are row or column vectors
    # and perform the dot product for each unit vector pair
    # ndim = 6
    # axis = xforms.check_axis(star_pair, ndim, axis=axis)

    if axis is 0:
        star1 = star_pair[0:3, :]
        star2 = star_pair[3:6, :]
        dotp = xforms.vector_dot(star1, star2)
    elif axis is 1:
        star1 = star_pair[:, 0:3]
        star2 = star_pair[:, 3:6]
        dotp = xforms.vector_dot(star1.T, star2.T)
    else:  # possibly an empty array
        raise ValueError("ERROR ["+str(__name__)+"]: {0} is not a valid axis. Error passed to {1}".format(axis, __name__))

    # sometimes the dot product returns values just over 1.0, so set to 1.0
    [print(x) for x in dotp if (x > 1.0005) ]
    dotp = np.array([x if x <= 1.0 else 1.0 for x in dotp])
    # return the arccosine to calculate the interstar angle
    return np.arccos(dotp)


def enhanced_pattern_shifting(candidateIdx):
    import numpy as np
    '''
    # DESCRIPTION: Creates patterned kernel for reducing the time spent from
    # iteration to iteration on the same points
    #
    # INPUTS: candidateIdx - ordered candidate star indices from list
    #
    #
    # OUTPUT: pKernel - patterned kernel for output into star ID algorithm
    #
    # Reference: Arnas, Fiahlo, and Mortari, "Fast and robust kernel generators
    # for star trackers" (2017).
    '''

    n = len(candidateIdx)
    if n < 3:
        print('Length less than minimum')
        return None
    u, c = np.unique(candidateIdx, return_counts=True)
    # check if duplicates exist
    dup = u[c > 1]
    if dup.size > 0:
        print('Duplicates found in candidate indices')
        return None
    pKernel = np.zeros([nchoosek(n, 3), 3], dtype=int)
    m = 0

    for dj in np.arange(1, n-2+1, dtype=int):
        for dk in np.arange(1, n-dj-1+1, dtype=int):
            for ii in [1, 2, 3]:
                for i in np.arange(ii, n-dj-dk+1, 3, dtype=int):
                    j = i + dj
                    k = j + dk
                    pKernel[m, :] = [i, j, k]
                    m = m + 1
    #          if m >= 20:
    #            return pKernel-1;
    return pKernel-1


def kvec_values(input_cat):
    # Calculates the k-vector slope and intercept of the
    # line given the input catalog
    import numpy as np
    if type(input_cat) is not np.ndarray:
        raise TypeError("ERROR ["+str(__name__)+"]: Input array is not ndarray in {}".format(__name__))

    input_cat = input_cat.flatten()
    s_idx = np.argsort(input_cat)
    sorted_cat = np.vstack((s_idx, input_cat[s_idx])).transpose()

    nrow = len(s_idx)
    ymin, ymax = (input_cat[s_idx[0]], input_cat[s_idx[-1]])
    # machine_prec = 2.220446049250313e-16 on my machine
    machine_prec = np.finfo(float).eps
    xi = machine_prec*max(abs(ymin), abs(ymax))
    m = (ymax-ymin+2*xi)/(nrow-1) if nrow > 1 else 0
    q = ymin-m-xi
    return m, q, sorted_cat


def ksearch(kvec, xin, d_thresh, d_cat, m, q):
    # Searches the catalog using interpolation technique
    # Reference: "k-vector range searching techniques", by D. Mortari and B. Neta
    import numpy as np
    ya = xin - d_thresh  # lower limit of search range
    yb = xin + d_thresh  # upper limit of search range

    ncol = len(d_cat[0])-1

    jbot = int(np.floor((ya-q)/m))
    jtop = int(np.ceil((yb-q)/m))

    # Guarantee that jbot and jtop remain within bounds of kvector
    jbot = min(max(jbot, 0), len(kvec)-1)
    jtop = min(jtop, len(kvec)-1)

    # Get the values of the angle between the two indices
    d = d_cat[int(kvec[jbot]):int(kvec[jtop]-1), :]
    # Remove the values that fall outside the range; see paper
    vfound = np.where(
        np.logical_and(
            d[:, ncol]+d_thresh >= xin,
            d[:, ncol]-d_thresh <= xin))[0]
    matchid = np.array(d_cat[vfound+int(kvec[jbot]), 0]).astype(int)
#    dist = d[vfound,1]-xin
    return matchid


def full_obs_match(x_obs, x_cat, isa_thresh):
    import numpy as np
    from array_transformations import vector_norm
    from math import sin, cos
    n_obs = len(x_obs[0])

    # initialize array of all entries -1: some IDs are 0
    idmatch = np.zeros([n_obs, 1],dtype=int)-1

    nmatches = 0
    cos_isa_thresh = cos(isa_thresh)

    for count, obs in enumerate(x_obs.transpose()):
        # calculate the distance between each vector in x_cat and the observation
        x_dot = np.dot(x_cat.T, obs)
        # find the vector(s) that fall under the threshold
        idx = np.where(x_dot >= cos_isa_thresh)[0]

        # original code below vvvvvv
        # if the number of vector distances below the threshold is one, add it to the number of matches
        #if idx.size == 1:
        #    idmatch[count] = idx
        #    nmatches += 1
        # original code above ^^^^^^^^^^

        # new code below vvvvvv
        #if there's only one return within the threshold above, then just use it
        if idx.size == 1:
            idmatch[count] = idx
            nmatches += 1

        #if there's more than one return, return the one with the smallest delta angle
        #TODO: FIXME: verify this is, indeed, grabbing the smallest delta angle as the above sure implies it's not?
        elif idx.size > 1:
            best_match = idx[0]
            for the_idx in idx:
                if x_dot[the_idx] < x_dot[best_match]: best_match=the_idx
            idmatch[count] = best_match
            nmatches += 1
        # new code above ^^^^^^^^^^

    # idmatch: mapping of stars from x_obs to index of x_cat
    # nmatches: total number of matched stars
    return idmatch, nmatches


def attitude_svd(ei, es):
    import numpy as np
    """
    Solves Wahba's problem using the SVD method.

    Parameters
    ----------
    ei :
        Catalog unit vectors in the inertial frame, 3xn matrix
    es :
        Measured unit vectors in the sensor frame, 3xn matrix

    Returns
    -------
    T :
        Rotation matrix from inertial frame to sensor frame

    G. Wahba, "A Least Squares Estimate of Spacecraft Attitude", (1965)
    """
    if type(ei) is not np.ndarray:
        raise TypeError("ERROR ["+str(__name__)+"]: Input 'ei' is not ndarray in {}".format(__name__))
    if type(es) is not np.ndarray:
        raise TypeError("ERROR ["+str(__name__)+"]: Input 'es' is not ndarray in {}".format(__name__))
    # Resize the es array
    # size = es.shape[1]
    # Transpose ei
    ei = np.transpose(ei)
    # Create array B
    B = np.zeros([3, 3], dtype=int)
    for index, row in enumerate(ei):
        # Columns of es
        # col = es[:, index].reshape((1, size))
        col = es[:, index]
        # Columns of es times rows of ei
        B = B + np.outer(col, row)
    # Perform SVD
    U, s, V = np.linalg.svd(B)
    d = np.linalg.det(U)*np.linalg.det(V)
    M = np.array([[1, 0, 0], [0, 1, 0], [0, 0, d]])
    T = U.dot(M).dot(V)
    return T


@timing_decorator
def triangle_isa_id(x_obs, x_cat, idx_star_pairs, isa_thresh, nmatch,
                    k, k_vector_interp, watchdog=None, verbose=False):
    import time
    import numpy as np
    import array_transformations as xforms

    # initialize values
    start_time = time.time()
    if watchdog is None: watchdog = 3600 #if no watchdog was provided, set it to an hour
    nmatches = 0
    nmatch_array = [] #used for troubleshooting
    unsolved_potentials = 0 #used for troubleshooting
    solved_potentials = 0 #used for troubleshooting
    idmatch = np.zeros(len(x_obs[0]))
    q_est = np.empty((4, ))
    q_est.fill(np.nan)
    # unpack star pairs and interstar angles from catalog
    star_pairs = idx_star_pairs[:, 0:2].astype(int)
    isa_cat = idx_star_pairs[:, 2:]
    # unpack k-vector interpolation values (slope and intercept)
    m, q = (k_vector_interp[0], k_vector_interp[1])

    # perform enhanced pattern shifting to produce star triplets
    n_obs = len(x_obs[0])
    if n_obs < nmatch:
        if verbose: print('Insufficient number of candidates to process for attitude estimation')
        return q_est, idmatch, nmatches
    mm = enhanced_pattern_shifting(np.arange(0, n_obs, 1, dtype=int))

    if verbose: print("    [triangle_isa_id]: length of enhanced pattern shifting return: " + str(len(mm)))
    for pair in np.arange(0, len(mm)-1, 1, dtype=int):
        t_idx = mm[pair]
        # calculate the interstar angle between each star pair with in the triplet
        star_pair_ab = np.concatenate(([x_obs[:, t_idx[0]], x_obs[:, t_idx[1]]])).reshape(1, 6)
        star_pair_bc = np.concatenate(([x_obs[:, t_idx[1]], x_obs[:, t_idx[2]]])).reshape(1, 6)
        star_pair_ac = np.concatenate(([x_obs[:, t_idx[0]], x_obs[:, t_idx[2]]])).reshape(1, 6)
        isa_ab = interstar_angle(star_pair_ab)
        isa_bc = interstar_angle(star_pair_bc)
        isa_ac = interstar_angle(star_pair_ac)

        # Find the indices & angles of the catalog entries within a
        # specified tolerance
        idx_ab = ksearch(k, isa_ab, isa_thresh, isa_cat, m, q)
        idx_bc = ksearch(k, isa_bc, isa_thresh, isa_cat, m, q)
        idx_ac = ksearch(k, isa_ac, isa_thresh, isa_cat, m, q)

        # Get the pair IDs from the catalog for each leg,
        # and repeat the match lists so that pair order matters
        pairs_match_ab = np.concatenate((star_pairs[idx_ab, :], star_pairs[idx_ab, ::-1]))
        pairs_match_bc = np.concatenate((star_pairs[idx_bc, :], star_pairs[idx_bc, ::-1]))
        pairs_match_ac = np.concatenate((star_pairs[idx_ac, :], star_pairs[idx_ac, ::-1]))

        for ii_ab in pairs_match_ab:
            A = ii_ab[0]
            B = ii_ab[1]

            # Find possible values for C from leg AC
            match_row_ac = np.where(pairs_match_ac[:, 0] == A)[0]
            # Find possible values for C from leg BC
            match_row_bc = np.where(pairs_match_bc[:, 0] == B)[0]

            if (time.time()-start_time) > watchdog:
                nmatches = 0
                idmatch = np.zeros(len(x_obs[0]))
                q_est = np.empty((4, ))
                q_est.fill(np.nan)
                if verbose: print("    [triangle_isa_id]: hit watchdog timeout (" + str(time.time()-start_time)+"s of "+str(watchdog)+"s")
                return q_est, idmatch, nmatches

            # if we find that the number of matches above are both 1,
            # check if those two corresponding C values are the same
            #TODO: this explicitly requires only 1 match, look into making this more robust and/or a way to find all 3 catalog angle matches in one go
            if match_row_ac.__len__() > 1 or match_row_bc.__len__() > 1:
                if match_row_ac.__len__() > 0 and match_row_bc.__len__() > 0:
                    #print(match_row_ac.__len__())
                    #print(match_row_bc.__len__())
                    #print('----')
                    unsolved_potentials+=1

            if match_row_ac.__len__() == 1 and match_row_bc.__len__() == 1:
                solved_potentials+=1


            #TODO: could ditch the above, BUT would then need to clean up lists by removing those that aren't common to all
                C1 = pairs_match_ac[match_row_ac, 1]
                C2 = pairs_match_bc[match_row_bc, 1]

                if C1 == C2:
                    # get
                    p = x_cat[:, [A, B, C1[0]]]
                    y = x_obs[:, mm[pair, :]]

                    # perform singular value decomposition using Wahba's problem
                    t_hat = attitude_svd(p, y)

                    # rotate all stars in the camera frame into the estimated attitude frame
                    x_obs_CATFRAME = xforms.vector_array_transform(t_hat.T, x_obs)

                    # check the number of stars that align with candidate stars
                    idmatch, nmatches = full_obs_match(x_obs_CATFRAME, x_cat, isa_thresh)

                    if verbose: nmatch_array+=[nmatches]

                    # if we meet or exceed the number of matches required, return the result
                    if nmatches >= nmatch:

                        if verbose:
                            print("    [triangle_isa_id]: Success-- found "+str(nmatches)+" matches")
                            print("    [triangle_isa_id]: nmatch array len: " +str(len(nmatch_array)))
                            print("    [triangle_isa_id]: number of nmatch=0: " +str(nmatch_array.count(0)))
                            print("    [triangle_isa_id]: number of nmatch=1: " +str(nmatch_array.count(1)))
                            print("    [triangle_isa_id]: number of nmatch=2: " +str(nmatch_array.count(2)))
                            print("    [triangle_isa_id]: number of nmatch=3: " +str(nmatch_array.count(3)))
                            print("    [triangle_isa_id]: number of nmatch=4: " +str(nmatch_array.count(4)))
                            print("    [triangle_isa_id]: number of nmatch=5: " +str(nmatch_array.count(5)))
                            print("    [triangle_isa_id]: unsolved_potentials: " +str(unsolved_potentials))
                            print("    [triangle_isa_id]: solved_potentials: " +str(solved_potentials))


                        # convert attitude matrix to quaternion
                        q_est = xforms.attitude_matrix2quat(t_hat)

                        return q_est, idmatch, nmatches

    # if no attitude found, return an attitude matrix of nan's
    if verbose:
        print("    [triangle_isa_id]: Failed-- found "+str(nmatches)+" matches")
        print("    [triangle_isa_id]: nmatch array len: " +str(len(nmatch_array)))
        print("    [triangle_isa_id]: number of nmatch=0: " +str(nmatch_array.count(0)))
        print("    [triangle_isa_id]: number of nmatch=1: " +str(nmatch_array.count(1)))
        print("    [triangle_isa_id]: number of nmatch=2: " +str(nmatch_array.count(2)))
        print("    [triangle_isa_id]: number of nmatch=3: " +str(nmatch_array.count(3)))
        print("    [triangle_isa_id]: number of nmatch=4: " +str(nmatch_array.count(4)))
        print("    [triangle_isa_id]: number of nmatch=5: " +str(nmatch_array.count(5)))
        print("    [triangle_isa_id]: unsolved_potentials: " +str(unsolved_potentials))
        print("    [triangle_isa_id]: solved_potentials: " +str(solved_potentials))

    nmatches = 0
    idmatch = np.zeros(len(x_obs[0]))
    q_est = np.empty((4, ))
    q_est.fill(np.nan)

    return q_est, idmatch, nmatches


def calculate_center_intensity(img, stats, min_star_area, max_star_area):
    import cv2 as cv
    import numpy as np

    nrows, ncols = img.shape

    # Find areas greater than minimum threshold and less than maximum
    min_idx = np.where(stats[:, cv.CC_STAT_AREA] >= min_star_area)[0]
    max_idx = np.where(stats[min_idx, cv.CC_STAT_AREA] <= max_star_area)[0]
    # Of those with at least the min threshold area,
    # find those whose area does not exceed max threshold
    idx = min_idx[max_idx]
    stats = stats[idx, :]

    # initialize center of intensity (COI)
    coi = np.zeros(shape=(len(stats), 2), dtype=float)
    intensities = np.zeros(shape=(len(stats), 1), dtype=float)

    #TODO: FIXME.  The below list switch is only due to the intensity sum issue detailed below
    coi = []

    for count, label in enumerate(stats):
        top_pixel    = label[cv.CC_STAT_TOP]
        left_pixel   = label[cv.CC_STAT_LEFT]

        # ensure that bottom and right pixels don't go outside image bounds
        bottom_pixel = min(nrows-1, label[cv.CC_STAT_HEIGHT] +  top_pixel)
        right_pixel  = min(ncols-1, label[cv.CC_STAT_WIDTH]  + left_pixel)
        # mask of intensities in ROI, in [row, col]
        intensity = img[top_pixel:bottom_pixel+1, left_pixel:right_pixel+1]
        intensity_sum = np.sum(intensity)
        intensities[count, :] = intensity_sum

        #TODO: FIXME: intensity sum can sometimes be 0, which kills everything.  How can it be 0?!

        #print('-----')
        #print(intensity)
        #print(intensity_sum)


        if intensity_sum > 0: #else is only here due to issue w/ intensity_count == 0

            # print(sum(arange(left_pixel, right_pixel+1)*sum(intensities, axis=0).T))
            x_bar = np.sum(np.arange(left_pixel, right_pixel+1)*np.sum(intensity, axis=0).T)/intensity_sum
            y_bar = np.sum(np.arange(top_pixel, bottom_pixel+1)*np.sum(intensity, axis=1).T)/intensity_sum
            #coi[count, :] = [x_bar, y_bar]
            coi+=[[x_bar, y_bar]]

    #these array concat methods are only here due to issue w/ intensity_count == 0
    intensities = intensities[intensities != 0]
    coi = np.array(coi)

    return coi, intensities


def undistort_image(img, idxUndistorted, idxDistorted, interpWeights):
    import array_transformations as xforms
    import numpy as np
    img_size = np.array(img.shape)
    # img = img.astype(float)/255.0
    img_u = np.zeros(img_size)
    # TODO: parallelize this code; it takes too long
    # ind2sub returns tuples
    # indexing = 0
    order = 'F'
    output_fmt = 'P'
    idxUndistorted_rc = xforms.ind2sub(
        img_size, idxUndistorted,order=order, output_fmt=output_fmt)
    idxDistorted_rc0 = xforms.ind2sub(
        img_size, idxDistorted[:, 0],order=order, output_fmt=output_fmt)
    idxDistorted_rc1 = xforms.ind2sub(
        img_size, idxDistorted[:, 1],order=order, output_fmt=output_fmt)
    idxDistorted_rc2 = xforms.ind2sub(
        img_size, idxDistorted[:, 2],order=order, output_fmt=output_fmt)
    idxDistorted_rc3 = xforms.ind2sub(
        img_size, idxDistorted[:, 3],order=order, output_fmt=output_fmt)

    w = interpWeights[:,0]*img[idxDistorted_rc0[0],idxDistorted_rc0[1]]
    x = interpWeights[:,1]*img[idxDistorted_rc1[0],idxDistorted_rc1[1]]
    y = interpWeights[:,2]*img[idxDistorted_rc2[0],idxDistorted_rc2[1]]
    z = interpWeights[:,3]*img[idxDistorted_rc3[0],idxDistorted_rc3[1]]
    img_u[idxUndistorted_rc[0],idxUndistorted_rc[1]] = (w + x + y + z)[:, np.newaxis]

    return round(img_u).astype(np.uint8)


def input_parser(caller_name):
    import sys

    argv = sys.argv
    if len(argv) > 1:
        argv = argv[1:]
    else:
        argv = []

    # When --help or no args are given, print this help
    usage_text = (
        "Run script using the following bash command:\n   "
        + caller_name + " -- [options]"
    )
    # print('**** {0} ****'.format(argv))
    return argv, usage_text
