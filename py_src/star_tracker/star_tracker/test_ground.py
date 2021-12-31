#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""

TBR

"""

################################
#LOAD LIBRARIES
################################
from unittest import TestCase
import ground
import numpy as np


class Test(TestCase):

    # def test_kvector(self):
    #     self.fail()
    #
    def test_equatorial2vector(self):
        reqd_dim = 2
        output_axis = 0
        n_output_cols_rows = 3
        nrows, ncols = (reqd_dim, 4)
        a = np.zeros((nrows, ncols))
        arr_check = np.zeros((n_output_cols_rows, ncols))
        arr_check[0:1, :] += np.ones((1, ncols))
        ra = np.zeros((4, ))
        de = np.zeros((4, ))

        self.assertIsNone(np.testing.assert_array_equal(
            ground.equatorial2vector(a, axis=None)[0], arr_check))
        self.assertIsNone(np.testing.assert_array_equal(
            ground.equatorial2vector(a, axis=None)[1], ra))
        self.assertIsNone(np.testing.assert_array_equal(
            ground.equatorial2vector(a)[1], ra))
        self.assertIsNone(np.testing.assert_array_equal(
            ground.equatorial2vector(a, axis=None)[2], de))
        with self.assertRaises(ValueError):
            ground.equatorial2vector(a, axis=1)
    #
    # def test_vector2equatorial(self):
    #     self.fail()
    #
    def test_lpq_orthonormal_basis(self):
        length = 3
        arr = np.zeros((length, 2))
        arr = np.array([[0, 0], [90, 0], [0, 90]])
        arr_check = np.zeros((3, ))

        los, p_hat, q_hat = ground.lpq_orthonormal_basis(arr)
        self.assertIsNone(
            np.testing.assert_array_equal(
                ground.lpq_orthonormal_basis(arr)[0], los))
    #
    # def test_proper_motion_correction(self):
    #     self.fail()
    #
    # def test_get_hipparcos_data(self):
    #     self.fail()
    #
    # def test_read_star_catalog(self):
    #     self.fail()
    #
    # def test_create_star_catalog(self):
    #     self.fail()
    #
    # def test_read_dir_images(self):
    #     self.fail()
    #
    # def test_create_darkframe(self):
    #     self.fail()
    #
    # def test_create_distortion_map(self):
    #     self.fail()
    #
    # def test_find_file(self):
    #     self.fail()
    #
    # def test_find_files_pattern(self):
    #     self.fail()
