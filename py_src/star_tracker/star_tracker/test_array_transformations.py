#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""

TBR

"""

################################
#LOAD LIBRARIES
################################
from unittest import TestCase
import array_transformations as xforms


class Test(TestCase):
    # def test_rx_passive(self):
    #     self.fail()
    #
    # def test_ry_passive(self):
    #     self.fail()
    #
    # def test_rz_passive(self):
    #     self.fail()
    #
    # def test_attitude_matrix2quat(self):
    #     self.fail()
    #
    # def test_quat2attitude_matrix(self):
    #     self.fail()
    #
    # def test_vector_dot(self):
    #     self.fail()
    #
    # def test_vector_norm(self):
    #     self.fail()
    #
    # def test_normalize_vector_array(self):
    #     self.fail()
    #
    def test_camera2homogeneous(self):
        self.fail()

    # def test_camera2vector(self):
    #     self.fail()
    #
    # def test_pixel2vector(self):
    #     self.fail()
    #
    # def test_vector2pixel(self):
    #     self.fail()
    #
    # def test_vector2homogeneous(self):
    #     self.fail()
    #
    # def test_homogeneous2vector(self):
    #     self.fail()
    #
    # def test_vector_array_transform(self):
    #     self.fail()
    #
    # def test_matrix_multiplication(self):
    #     self.fail()
    #
    # def test_sub_ind_format(self):
    #     self.fail()
    #
    # def test_sub2ind(self):
    #     self.fail()
    #
    # def test_ind2sub(self):
    #     self.fail()

    def test_check_axis(self):
        import numpy as np
        reqd_dim = 2
        output_axis = 1
        nrows, ncols = (reqd_dim, 1)
        a = np.arange(nrows*ncols)
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=None), output_axis)
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=output_axis), output_axis)
        with self.assertRaises(ValueError):
            xforms.check_axis(a, reqd_dim, axis=0)

        nrows, ncols = (1, reqd_dim)
        output_axis = 1
        a = np.arange(nrows*ncols)
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=None), output_axis)
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=output_axis), output_axis)
        with self.assertRaises(ValueError):
            xforms.check_axis(a, reqd_dim, axis=0)

        reqd_dim = 3
        nrows, ncols = (reqd_dim, 4)
        output_axis = 0
        a = np.arange(nrows*ncols).reshape((nrows, ncols))
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=None), output_axis)
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=output_axis), output_axis)
        with self.assertRaises(ValueError):
            xforms.check_axis(a, reqd_dim, axis=1)

        nrows, ncols = (4, reqd_dim)
        output_axis = 1
        a = np.arange(nrows*ncols).reshape((nrows, ncols))
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=None), output_axis)
        self.assertEqual(
            xforms.check_axis(a, reqd_dim, axis=output_axis), output_axis)
        with self.assertRaises(ValueError):
            xforms.check_axis(a, reqd_dim, axis=0)

        nrows, ncols = (0, 0)
        a = np.arange(nrows*ncols).reshape((nrows, ncols))
        with self.assertRaises(ValueError):
            xforms.check_axis(a, reqd_dim, axis=None)
            xforms.check_axis(a, reqd_dim, axis=0)
            xforms.check_axis(a, reqd_dim, axis=1)

        with self.assertRaises(ValueError):
            # dimensions exceed 2
            a = np.arange(nrows*ncols).reshape((nrows, ncols, 1))
            xforms.check_axis(a, reqd_dim, axis=None)
            xforms.check_axis(a, 0, axis=None)

        with self.assertRaises(TypeError):
            # not supplied with numpy array as arr
            xforms.check_axis((0, 0), reqd_dim, axis=None)

################################
#MAIN CODE
################################
if __name__ == '__main__':
    import unittest
    unittest.main()
