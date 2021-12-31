#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
array_transformations.py

File containing functions used for matrix transformations

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
from numpy import pi

rad2deg = 180/pi
deg2rad = pi/180


################################
#SUPPORT FUNCTIONS
################################
def check_axis_decorator(dim, verbose=False):
    def outer_decorator(func):
        def axis_check(*args, **kwargs):
            axis = kwargs['axis'] if 'axis' in kwargs.keys() else None
            if verbose:
                print("Chosen axis is {}".format(axis))
            kwargs['axis'] = check_axis(args[0], dim, axis)
            return func(*args, **kwargs)
        return axis_check
    return outer_decorator


def rx_passive(theta):
    import numpy as np
    r = np.array([
        [1, 0, 0],
        [0, np.cos(theta), np.sin(theta)],
        [0, -np.sin(theta), np.cos(theta)]])
    return r


def ry_passive(theta):
    import numpy as np
    r = np.array([
        [np.cos(theta), 0, -np.sin(theta)],
        [0, 1, 0],
        [np.sin(theta), 0, np.cos(theta)]])
    return r


def rz_passive(theta):
    import numpy as np
    r = np.array([
        [np.cos(theta), np.sin(theta), 0],
        [-np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]])
    return r


def attitude_matrix2quat(a):
    from numpy import trace, argmax, array, linalg
    tra = trace(a)
    a11 = a[0, 0]
    a22 = a[1, 1]
    a33 = a[2, 2]

    idx = argmax([abs(a11), abs(a22), abs(a33), abs(tra)])

    q = {
        0: array([1 + 2*a11 - tra,
                  a[0, 1] + a[1, 0],
                  a[0, 2] + a[2, 0],
                  a[1, 2] - a[2, 1]]),
        1: array([a[1, 0] + a[0, 1],
                  1 + 2*a22 - tra,
                  a[1, 2] + a[2, 1],
                  a[2, 0] - a[0, 2]]),
        2: array([a[2, 0] + a[0, 2],
                  a[2, 1] + a[1, 2],
                  1 + 2*a33 - tra,
                  a[0, 1] - a[1, 0]]),
        3: array([a[1, 2] - a[2, 1],
                  a[2, 0] - a[0, 2],
                  a[0, 1] - a[1, 0],
                  1 + tra])
        # nan is default if idx not found
    }.get(int(idx), [float('nan'), float('nan'), float('nan'), float('nan')])
    q = q/linalg.norm(q)
    return q


def quat2attitude_matrix(q13, q4):
    import numpy as np
    # Calculate norm of the quaternion
    q1 = q13[0]
    q2 = q13[1]
    q3 = q13[2]
    qnorm2 = q1**2 + q2**2 + q3**2 + q4**2
    # Cross product matrix
    # See Eq 2.55 from [Markley & Crassidis, 2014]
    # q13x = [   0, -q3,  q2
    # 	        q3,   0, -q1
    # 		   -q2,  q1,   0  ]

    # Quaternion representation
    # See Eq 2.125 from [Markley & Crassidis, 2014]
    # a = (q4^2-norm(q13)^2)*eye(3) - 2*q4*q13x + 2*q13*q13';
    a = np.array([[q1**2-q2**2-q3**2+q4**2, 2*(q1*q2+q3*q4),         2*(q1*q3-q2*q4)],
                  [2*(q2*q1-q3*q4),        -q1**2+q2**2-q3**2+q4**2, 2*(q2*q3+q1*q4)],
                  [2*(q3*q1+q2*q4),         2*(q3*q2-q1*q4),        -q1**2-q2**2+q3**2+q4**2]])
    a = a/qnorm2
    return a


def vector_dot(a, b):
    # Calculates the dot product of two arrays of vectors using einstein summation
    # Expects: a-> nxn, b-> nxm
    import numpy as np
    assert len(a) == len(b)
    return np.einsum("ij..., ij...->j...", a, b)


def vector_norm(v):
    # Calculates the vector norm of an array of N vectors using einstein summation in vector_dot function
    # from math import sqrt
    import numpy as np
    if v.ndim <= 1:
        v = v.reshape(len(v), 1)
    return np.sqrt(vector_dot(v, v))


def normalize_vector_array(v):
    # Normalizes each vector in array of N vectors
    return v / vector_norm(v)


@check_axis_decorator(2)
def camera2homogeneous(p, axis=None):
    # Accepts (2, n) or (n, 2) arrays
    # Converts camera (x, y) coordinates to homogeneous space; returns (3,n) numpy array
    # this looks like (x, y, z) for z=1, which is what we use
    import numpy as np
    nvec = p.shape[1] if axis == 0 else p.shape[0]
    p_homogeneous = np.vstack((p, np.ones((1, nvec))))
    return p_homogeneous


def camera2vector(xy):
    # Converts camera (x, y, 1) coordinates to unit vectors
    # assert len(xy)
    return normalize_vector_array(camera2homogeneous(xy))


def pixel2vector(cinv, p):
    # Transforms pixel coordinates (u ,v) to (x, y) coordinates in camera frame
    # h = pixel2homogeneous(p)
    h = camera2homogeneous(p, axis=0)
    xy_unnormalized = vector_array_transform(cinv, h)
    return normalize_vector_array(xy_unnormalized)


def vector2pixel(c, v):
    # TODO: add checks/input verification
    return vector2homogeneous(c, v)[0:2, :]


def vector2homogeneous(c, v):
    # TODO: add checks/input verification
    return vector_array_transform(c, v)/v[2, :]  # [0:2, :]


def homogeneous2vector(cinv, h):
    # TODO: add checks/input verification
    # Transforms homogeneous coordinates (u ,v, 1) to (x, y) coordinates in camera frame
    return normalize_vector_array(vector_array_transform(cinv, h))


def vector_array_transform(transform_matrix, v):
    # TODO: add checks/input verification
    # Performs matrix multiplication using numpy arrays rather than using matrix
    import numpy as np
    if len(v.shape) is 1:
        v = v[np.newaxis].T
    # assert len(transform_matrix[0]) == len(v)
    return np.einsum('ij,jk->ik', transform_matrix, v)


def matrix_multiplication(mat1, *argv):
    mat_result = mat1
    if argv is None:
        return mat_result
    for mat in argv:
        if len(mat_result[0]) == len(mat):
            raise ValueError("ERROR ["+str(__name__)+"]: Array sizes for operands are incompatible")
        mat_result = vector_array_transform(mat_result, mat)
    return mat_result


def sub_ind_format(output_fmt):
    if output_fmt is not 'M' and output_fmt is not 'P':
        raise ValueError("ERROR ["+str(__name__)+"]: Invalid input for format type")
    return 0 if output_fmt is 'P' else 1



def sub2ind(img_size, row, col, order='F', indexing=0, output_fmt='P'):
    import numpy as np
    # img_size: (rows,cols) tuple of image size
    # row: (n,) array of column indices
    # col: (n,) array of row indices
    # returns (n,) array of indices
    # nrows = img_size[0]
    # ncols = img_size[1]
    output_idx = sub_ind_format(output_fmt)
    rc = np.array([row, col])-indexing
    return np.ravel_multi_index(rc, np.flip(img_size), order=order) + output_idx
    # return array((col-1) * nrows + row, dtype=int)


def ind2sub(img_size, idx, order='F', output_fmt='P'):
    # img_size: (rows,cols) tuple of image size
    # idx: (n,) array of image indices
    # returns (2,n) array of row and column indices
    import numpy as np
    # nrows = img_size[0]
    # # ncols = img_size[1]
    # col = idx // nrows + 1
    # row = idx % nrows + 1
    output_idx = np.array(sub_ind_format(output_fmt))
    return np.unravel_index(idx, img_size, order=order) + output_idx
    # return array([row, col], dtype=int)


def check_axis(arr, dim, axis=None):
    # ex: if arr is a 2x3 array, and the array has pairs of numbers (x,y)
    # so the dim=2
    # 2x3 => axis = 0
    #  | x1 x2 x3 |
    #  | y1 y2 y3 |
    import numpy as np
    if type(arr) is not np.ndarray:
        raise TypeError("ERROR ["+str(__name__)+"]: Input array is not ndarray")
    # determine/verify if star pairs are row or column vectors according to axis
    shape = arr.shape
    ndims = arr.ndim
    # arr could be a 1D array, and
    nrows, ncols = (0, shape[0]) if ndims == 1 else (shape[0], shape[1])
    if ndims > 2:
        raise ValueError("ERROR ["+str(__name__)+"]: Input array has more than two dimensions "+str(ndims))
    if shape[0] is 0:
        raise ValueError("ERROR ["+str(__name__)+"]: Numpy array is empty")
    if (nrows != dim and ncols != dim) or dim < 1:
        raise ValueError("ERROR ["+str(__name__)+"]: Length of vector array must by nx2 or 2xn")
    if axis is None:
        if nrows == ncols:
            raise ValueError("ERROR ["+str(__name__)+"]: Length is equal for star pairs and no axis provided")
        # TODO: verify that every function returns the same axis \\
        #  according to this convention
        axis = 0 if nrows == dim else 1
    if axis == 0 and nrows != dim:
        raise ValueError("ERROR ["+str(__name__)+"]: Axis {}: Columns: {}, Reqd: {}".format(axis, ncols, dim))
    if axis == 1 and ncols != dim:
        raise ValueError("ERROR ["+str(__name__)+"]: Axis {}: Columns: {}, Reqd: {}".format(axis, nrows, dim))
    # axis = 0 if ncols == dim else 1
    # axis will be 0 if each column contains a vector
    # axis will be 1 if each row contains a vector
    return axis
