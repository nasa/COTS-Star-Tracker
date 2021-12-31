#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
cam_matrix.py

File containing functions used for camera calibration and transformations

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
#SUPPORT FUNCTIONS
################################
def cam_matrix_from_params(dx, dy, up, vp, sk):
    import numpy as np
    return np.array([[dx, sk, up], [0, dy, vp], [0, 0, 1]])


def cam_matrix_inv(c):
    import numpy as np
    dx = c[0, 0]
    dy = c[1, 1]
    sk = c[0, 1]
    up = c[0, 2]
    vp = c[1, 2]

    return np.array([
        [1/dx, -sk/(dx*dy), (sk*vp-dy*up)/(dx*dy)],
        [0, 1/dy, -vp/dy],
        [0, 0, 1]])


def cam2fov(cinv, nrow, ncol):
    import numpy as np
    import array_transformations as xforms

    # get (1,1) pixel in unit vector form
    p1 = xforms.vector_array_transform(cinv, np.ones((3, 1)))
    p1 = p1/xforms.vector_norm(p1)
    # get boresight unit vector
    p2 = xforms.vector_array_transform(cinv, np.array([[nrow/2], [ncol/2], [1]]))
    p2 = p2/xforms.vector_norm(p2)
    # take the arccos of the dot product of the two
    return np.arccos(xforms.vector_dot(p1, p2))


def read_cam_json(cam_config_file):
    import numpy as np
    import os
    import json
    print("")
    if os.path.exists(cam_config_file):
        # print("reading camera parameters from file...")
        with open(cam_config_file, "r") as read_file:

            data = json.load(read_file)

        cam_resolution = np.array([data["resolution"][0], data["resolution"][1]], dtype=int)  # pixels, [cols, rows]

        # TODO: verify that every coefficient has a corresponding value
        coeff_names = ['k1', 'k2', 'p1', 'p2', 'k3']
        k_1, k_2, p_1, p_2, k_3 = [data.get(x) for x in coeff_names]

        # we use fx/dx and fy/dy interchangeably
        try:
            dx = data["dx"]
            dy = data["dy"]
        except:
            dx = data["fx"]
            dy = data["fy"]

        up = data["up"]
        vp = data["vp"]
        skew = data["skew"]

    else: #TODO one day make assumptions based on input image
        print(cam_config_file)
        print('ERROR ['+str(__name__)+']: camera config file not found/failed to load!  Please fix the file and/or path and try again.  Exiting...')

        #pixel_pitch = [0.00112, 0.00112]  # mm
        #focal_length = 3  # mm
        #prncple_pt = [cam_resolution[0]/2, cam_resolution[1]/2]  # pixels
        #k_1, k_2, p_1, p_2, k_3 = (None, None, None, None, None)
        ## convert to opencv units
        #prncple_pt_cv2 = [prncple_pt[0]*pixel_pitch[0], prncple_pt[1]*pixel_pitch[1]]  # mm
        #focal_length_cv2 = [focal_length/pixel_pitch[0], focal_length/pixel_pitch[1]]  # pixels
        #dx = focal_length_cv2[0]
        #dy = focal_length_cv2[1]
        #up = prncple_pt_cv2[0]
        #vp = prncple_pt_cv2[1]
        #skew = 0

    # distortion coeffs
    dist_coefs = np.array([k_1, k_2, p_1, p_2, k_3])
    dist_coefs[dist_coefs is None] = 0
    # populate camera matrix
    camera_matrix = cam_matrix_from_params(dx, dy, up, vp, skew)

    return camera_matrix, cam_resolution, dist_coefs


def fov2sensorArray(fov, f, cam_resolution):
    #     field-of-view    [in deg ]
    #  ---no of pixels--   [integer]
    #           ---dx---   [unitless] = focal_length/pixel_pitch
    #  \       |       /
    #   \      |      /
    #    \     |theta/  [angle]
    #     \  f |    /   [in mm]
    #      \   |   /
    #       \  |  /
    #        \ | /
    #         \|/
    #          *
    #       pinhole
    import numpy as np
    # calculate the pixel sensor array size and pixel size
    sensor_size = 2*f*np.tan(np.radians(fov)/2)
    pixel_pitch = sensor_size/cam_resolution
    return sensor_size, pixel_pitch


def focallength2fov(f, sensor_size):
    import numpy as np
    return np.degrees(2*np.arctan(sensor_size/(2*f)))


def cam_matrix2sensorArray(camera_matrix, f, cam_resolution):
    import numpy as np
    #     field-of-view    [in deg ]
    #  ---no of pixels--   [integer]
    #           ---dx---   [unitless] = focal_length/pixel_pitch
    #  \       |       /
    #   \      |      /
    #    \     |theta/  [angle]
    #     \  f |    /   [in mm]
    #      \   |   /
    #       \  |  /
    #        \ | /
    #         \|/
    #          *
    #       pinhole
    # calculate the pixel sensor array size and pixel size
    focal_lengths = np.array([camera_matrix[0, 0], camera_matrix[1, 1]])
    pixel_pitch = f/focal_lengths
    sensor_size = pixel_pitch*np.array(cam_resolution).astype(float)
    fov = focallength2fov(f, sensor_size)
    return sensor_size, pixel_pitch, fov


