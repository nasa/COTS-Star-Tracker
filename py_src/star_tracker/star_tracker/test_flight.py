#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""

TBR

"""

################################
#LOAD LIBRARIES
################################
from unittest import TestCase
import flight
import numpy as np


class TestFlight(TestCase):
    def test_nchoosek(self):
        self.assertEqual(flight.nchoosek(1, 1), 1)
        self.assertEqual(flight.nchoosek(5, 2), 10)
        self.assertEqual(flight.nchoosek(5, 1), 5)
        self.assertEqual(flight.nchoosek(5, 10), 0)
        self.assertEqual(flight.nchoosek(50, 5), 2118760)
        self.assertEqual(flight.nchoosek(50, 10), 10272278170)
        self.assertEqual(flight.nchoosek(0, 1), 0)
        self.assertEqual(flight.nchoosek(0, 0), 1)
        self.assertEqual(flight.nchoosek(1, 0), 1)
        self.assertEqual(flight.nchoosek(-1, 1), 0)
        self.assertEqual(flight.nchoosek(1, -10), 0)
        # self.assertIsNone(flight.nchoosek(0, 1))
        # self.assertIsNone(flight.nchoosek(1, 0))
        # self.assertIsNone(flight.nchoosek(-1, 1))
        # self.assertIsNone(flight.nchoosek(1, -10))
    #
    def test_interstar_angle(self):
        arr = np.array([[1, 0, 0, 0, 1, 0]])
        # flight.interstar_angle(arr, axis=None)
        self.assertEqual(
            flight.interstar_angle(arr, axis=1), np.pi/2)
        arr = np.array([[1], [0], [0], [1], [0], [0]])
        self.assertEqual(
            flight.interstar_angle(arr, axis=0), 0)
        # add tests for corner cases and random vectors and vector arrays

        with self.assertRaises(ValueError):
            arr = np.array([[1, 0, 0, 0, 0]])
            flight.interstar_angle(arr, axis=None)
            arr = np.array([[1], [0], [0], [0], [0], [0], [0]])
            flight.interstar_angle(arr, axis=None)
            arr = np.array([[1], [0], [0], [0], [0], [0]])
            flight.interstar_angle(arr, axis=2)

    def test_enhanced_pattern_shifting(self):
        arr = np.arange(0, 4, 1)
        arr_check = np.array([[0, 1, 2],
                              [1, 2, 3],
                              [0, 1, 3],
                              [0, 2, 3]])
        self.assertIsNone(
            np.testing.assert_array_equal(
                flight.enhanced_pattern_shifting(arr), arr_check))
        arr = np.arange(0, 3, 1)
        arr_check = np.array([[0, 1, 2]])
        self.assertIsNone(
            np.testing.assert_array_equal(
                flight.enhanced_pattern_shifting(arr), arr_check))
        arr = np.arange(0, 6, 1)
        arr_check = np.array([[0, 1, 2],
                              [3, 4, 5],
                              [1, 2, 3],
                              [2, 3, 4],
                              [0, 1, 3],
                              [1, 2, 4],
                              [2, 3, 5],
                              [0, 1, 4],
                              [1, 2, 5],
                              [0, 1, 5],
                              [0, 2, 3],
                              [1, 3, 4],
                              [2, 4, 5],
                              [0, 2, 4],
                              [1, 3, 5],
                              [0, 2, 5],
                              [0, 3, 4],
                              [1, 4, 5],
                              [0, 3, 5],
                              [0, 4, 5]])
        self.assertIsNone(
            np.testing.assert_array_equal(
                flight.enhanced_pattern_shifting(arr), arr_check))

        arr = np.array([])
        self.assertIsNone(flight.enhanced_pattern_shifting(arr))
        arr = np.array([0])
        self.assertIsNone(flight.enhanced_pattern_shifting(arr))
        arr = np.array([10, 11])
        self.assertIsNone(flight.enhanced_pattern_shifting(arr))
        arr = np.array([10, 11, 12, 11, 15, 20])
        self.assertIsNone(flight.enhanced_pattern_shifting(arr))

    def test_kvec_values(self):
        # import ground
        arr = np.arange(0, 1000, 100).astype(float)/10.0
        # m, q = flight.kvec_values(istar_angle, istar_idx)
        # print("{}\n{}".format(m, q))
        m_test, q_test = (10.000000000000005, -10.000000000000025)
        self.assertAlmostEqual(flight.kvec_values(arr)[0], m_test)
        self.assertAlmostEqual(flight.kvec_values(arr)[1], q_test)
        arr = np.arange(0.0, 6.0, 1)
        m_test, q_test = (1.0000000000000004, -1.0000000000000016)
        self.assertAlmostEqual(flight.kvec_values(arr)[0], m_test)
        self.assertAlmostEqual(flight.kvec_values(arr)[1], q_test)
        arr = np.linspace(0.3, 0.6, num=6)
        arr[0] = 0.3247375897416183
        arr[-1] = 0.6320437092965802
        m_test, q_test = (0.06146122391099242, 0.26327636583062575)
        self.assertAlmostEqual(flight.kvec_values(arr)[0], m_test)
        self.assertAlmostEqual(flight.kvec_values(arr)[1], q_test)

    # def test_ksearch(self):
    #     self.fail()
    #
    def test_full_obs_match(self):
        # TODO: check for input types and incorrect data handling
        x_obs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        x_cat = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        isa_thresh = 0.0001
        nmatch_test = 2
        idmatch_test = np.array([[-1], [0], [2]])
        self.assertIsNone(
            np.testing.assert_array_equal(
                flight.full_obs_match(
                    x_obs, x_cat, isa_thresh)[0], idmatch_test))
        self.assertIsNone(
            np.testing.assert_array_equal(
                flight.full_obs_match(
                    x_obs, x_cat, isa_thresh)[1], nmatch_test))
    #
    def test_attitude_svd(self):
        # TODO: check for input types and incorrect data handling
        x_obs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        x_cat = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        t_test = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        # t = flight.attitude_svd(x_obs, x_cat)
        self.assertIsNone(
            np.testing.assert_array_equal(
                flight.attitude_svd(x_obs, x_cat), t_test))
    #
    # def test_triangle_isa_id(self):
    #     self.fail()
    #
    # def test_calculate_center_intensity(self):
    #     self.fail()
    #     img = np.zeros((128, 256))
    #     stats = None
    #     min_star_area = 5
    #     max_star_area = 30
    #     coi, intensities = flight.calculate_center_intensity(
    #         img, stats, min_star_area, max_star_area)
    #
    # def test_undistort_image(self):
    #     self.fail()
    #
    # def test_input_parser(self):
    #     self.fail()


if __name__ == '__main__':
    import unittest
    unittest.main()
