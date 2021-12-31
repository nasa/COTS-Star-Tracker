#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""

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
import flight
import numpy as np
import cam_matrix as cm
import array_transformations
from itertools import combinations

################################
#MAIN CODE
################################
# Load data
u = np.load('u.npy')
d_cat = np.load('d_cat.npy')
star_pairs = np.load('star_pairs.npy')
# Image
im_filename = "../data/UpSouth(6-29-19,0030-0040)/Image-1.bmp"
# Intrinsic camera parameters matrix
c = np.transpose(np.array([[7.891092810360696e+03, 0, 0], [0,
                           7.906007692160952e+03, 0], [1.274779317743486e+03,
                                                       1.066692381188306e+03,
                                                       1]]))
# Inverse camera matrix
cinv = cm.cam_matrix_inv(c)
# 2xn pixel point matrix
pixpoints = np.transpose(np.array([[1194, 822], [1323, 956], [1114, 933],
                                   [1502, 752], [962, 679], [1396, 1025],
                                   [1602, 957], [1026, 815], [1342, 682],
                                   [1563, 633]]))

# Unit vector of points in camera frame
vector = array_transformations.pixel2vector(cinv, pixpoints)
# Create array of pairs
N = len(vector[0])
c_star_pairs = np.array(list(combinations(np.arange(0, N, 1, dtype=int), 2)),
                        dtype=int)
# Create nx6 matrix of candidate star pairs
cand_sp = np.concatenate([vector[:, c_star_pairs[:, 0]],
                          vector[:, c_star_pairs[:, 1]]], axis=0)
# Find ISA
c_ISA = flight.interstar_angle(cand_sp)

# Now find stars
x_obs = vector
x_cat = u
d_thresh = 0.03
nmatch = 5

k = np.load('k.npy')
m = np.load('m.npy')
q = np.load('q.npy')
t_hat, idmatch, nmatches = flight.triangle_isa_id(x_obs, d_thresh, nmatch,
                                                  x_cat, d_cat, star_pairs, k,
                                                  m, q)

