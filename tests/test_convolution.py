#!/usr/bin/python3

from fractions import Fraction

import unittest

import nrconv

class TestEdgeCase(unittest.TestCase):
    def test_non_rectangular_convolution_edge_diagonal1(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)), (Fraction(7, 1), Fraction(7, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_edge(
            list1, list2, geometry, prime)
        want = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_edge_diagonal2(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(7, 1)),
                    (Fraction(7, 1), Fraction(0, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_edge(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_edge_halfdiagonal1(self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(7, 1), Fraction(7, 2))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_edge(
            list1, list2, geometry, prime)
        want = [1, 0, 0, 3, 0, 0, 5, 0, 0, 7, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_edge_halfdiagonal2(self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(7, 2), Fraction(7, 1)),
                    (Fraction(0, 1), Fraction(0, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_edge(
            list1, list2, geometry, prime)
        want = [1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_edge_vertical(self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(2, 1), Fraction(0, 1)),
                    (Fraction(2, 1), Fraction(7, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_edge(
            list1, list2, geometry, prime)
        want = [3, 3, 3, 3, 3, 3, 3, 3]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_edge_offset(self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(1, 1), Fraction(7, 2)),
                    (Fraction(7, 1), Fraction(1, 2))]
        prime = nrconv.create_ntt_prime(list1, list2)
        _, result = nrconv.convolution.non_rectangular_convolution_edge(
            list1, list2, geometry, prime)
        want = 2
        self.assertEqual(result, want)


class TestRectangleCase(unittest.TestCase):
    def test_non_rectangular_convolution_rectangle_full(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)), (Fraction(7, 1), Fraction(7, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = [1, 2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_int_part(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(1, 1), Fraction(2, 1)), (Fraction(4, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = [1, 2, 3, 4, 4, 3, 2, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_frac_part(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(3, 4), Fraction(13, 2)),
                    (Fraction(13, 3), Fraction(5, 3))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = [1, 2, 3, 4, 4, 3, 2, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_int_dot(self):
        list1 = [1, 1, 3, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 7, 1, 1]
        geometry = [(Fraction(2, 1), Fraction(5, 1)), (Fraction(2, 1), Fraction(5, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = [21]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_intfrac_dot(self):
        list1 = [1, 1, 3, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 7, 1, 1]
        geometry = [(Fraction(6, 3), Fraction(5, 1)),
                    (Fraction(6, 3), Fraction(5, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = [21]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_frac_dot(self):
        list1 = [1, 1, 3, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 7, 1, 1]
        geometry = [(Fraction(3, 2), Fraction(17, 4)),
                    (Fraction(5, 2), Fraction(16, 3))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = [21]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_frac_empty1(self):
        list1 = [1, 1, 3, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 7, 1, 1]
        geometry = [(Fraction(4, 3), Fraction(17, 8)),
                    (Fraction(5, 3), Fraction(16, 3))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = []
        self.assertEqual(result, want)
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_frac_empty2(self):
        list1 = [1, 1, 3, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 7, 1, 1]
        geometry = [(Fraction(17, 8), Fraction(4, 3)),
                    (Fraction(16, 3), Fraction(5, 3))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = []
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_frac_part_offset(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(3, 4), Fraction(13, 2)),
                    (Fraction(13, 3), Fraction(5, 3))]
        prime = nrconv.create_ntt_prime(list1, list2)
        _, result = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = 3
        self.assertEqual(result, want)


class TestAxisAlignedTriangleCase(unittest.TestCase):
    def test_non_rectangular_convolution_triangle_axis_aligned_degenerated1(
            self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(4, 1), Fraction(4, 1)), (Fraction(4, 1), Fraction(7, 1)), (Fraction(4, 1), Fraction(4, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = [5, 5, 5, 5]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_degenerated2(
            self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(4, 3), Fraction(8, 2)),
                    (Fraction(19, 4), Fraction(4, 1)),
                    (Fraction(8, 6), Fraction(12, 3))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = [3, 4, 5]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_degenerated3(
            self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(4, 3), Fraction(5, 2)),
                    (Fraction(19, 4), Fraction(10, 4)),
                    (Fraction(8, 6), Fraction(15, 6))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = []
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_degenerated_offset(
            self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(4, 3), Fraction(8, 2)),
                    (Fraction(19, 4), Fraction(4, 1)),
                    (Fraction(8, 6), Fraction(12, 3))]
        prime = nrconv.create_ntt_prime(list1, list2)
        _, result = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = 6
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_small1(self):
        list1 = [0, 1, 2, 3, 4, 5, 6, 7]
        list2 = [0, 1, 2, 3, 4, 5, 6, 7]
        geometry = [(Fraction(39, 10), Fraction(32, 10)),
                    (Fraction(41, 10), Fraction(28, 10)),
                    (Fraction(41, 10), Fraction(32, 10))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = [12]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_small2(self):
        list1 = [0, 1, 2, 3, 4, 5, 6, 7]
        list2 = [0, 1, 2, 3, 4, 5, 6, 7]
        geometry = [(Fraction(39, 10), Fraction(32, 10)),
                    (Fraction(42, 10), Fraction(28, 10)),
                    (Fraction(42, 10), Fraction(32, 10))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = []
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_small3(self):
        list1 = [0, 1, 2, 3, 4, 5, 6, 7]
        list2 = [0, 1, 2, 3, 4, 5, 6, 7]
        geometry = [(Fraction(39, 10), Fraction(32, 10)),
                    (Fraction(39, 10), Fraction(28, 10)),
                    (Fraction(41, 10), Fraction(28, 10))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = [12]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_small4(self):
        list1 = [0, 1, 2, 3, 4, 5, 6, 7]
        list2 = [0, 1, 2, 3, 4, 5, 6, 7]
        geometry = [(Fraction(39, 10), Fraction(32, 10)),
                    (Fraction(39, 10), Fraction(28, 10)),
                    (Fraction(42, 10), Fraction(28, 10))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = [12]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_small_offset(
            self):
        list1 = [0, 1, 2, 3, 4, 5, 6, 7]
        list2 = [0, 1, 2, 3, 4, 5, 6, 7]
        geometry = [(Fraction(39, 10), Fraction(32, 10)),
                    (Fraction(39, 10), Fraction(28, 10)),
                    (Fraction(41, 10), Fraction(28, 10))]
        prime = nrconv.create_ntt_prime(list1, list2)
        _, result = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = 7
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_large(self):
        list1 = [0, 1, 2, 3, 4, 5, 6, 7]
        list2 = [0, 1, 2, 3, 4, 5, 6, 7]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(6, 1), Fraction(6, 1)),
                    (Fraction(0, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = [0, 0, 1, 2, 7, 10, 22, 28, 43, 38, 49, 30, 36]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_axis_aligned_big_offset(
            self):
        list1 = [0, 1, 2, 3, 4, 5, 6, 7]
        list2 = [0, 1, 2, 3, 4, 5, 6, 7]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(6, 1), Fraction(6, 1)),
                    (Fraction(0, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        _, result = nrconv.convolution.non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, prime)
        want = 0
        self.assertEqual(result, want)


class TestAxisArbitraryTriangleCase(unittest.TestCase):
    def test_non_rectangular_convolution_triangle_case_1_1_1(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(4, 1), Fraction(2, 1)),
                    (Fraction(6, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_1_1_2(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(2, 1), Fraction(4, 1)),
                    (Fraction(6, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 0, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_1_2_1(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(6, 1)),
                    (Fraction(4, 1), Fraction(4, 1)),
                    (Fraction(6, 1), Fraction(0, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 0, 0, 0, 7, 4, 1, 0, 0, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_1_2_2(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(6, 1)),
                    (Fraction(2, 1), Fraction(2, 1)),
                    (Fraction(6, 1), Fraction(0, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 0, 1, 4, 7, 0, 0, 0, 0, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_1_2_1(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(3, 1), Fraction(6, 1)),
                    (Fraction(6, 1), Fraction(0, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [1, 1, 2, 3, 3, 4, 5, 3, 2, 1, 0, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_1_2_2(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(6, 1)),
                    (Fraction(3, 1), Fraction(0, 1)),
                    (Fraction(6, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 1, 2, 3, 5, 4, 3, 3, 2, 1, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_1_1_1(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(6, 1), Fraction(3, 1)),
                    (Fraction(0, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [1, 1, 2, 3, 3, 4, 5, 3, 2, 1, 0, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_1_1_2(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(6, 1), Fraction(0, 1)),
                    (Fraction(0, 1), Fraction(3, 1)),
                    (Fraction(6, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 1, 2, 3, 5, 4, 3, 3, 2, 1, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_2_1(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(6, 1), Fraction(3, 1)),
                    (Fraction(3, 1), Fraction(6, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 0, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_2_2(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(3, 1), Fraction(6, 1)),
                    (Fraction(6, 1), Fraction(3, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 0, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_2_3(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(6, 1), Fraction(6, 1)),
                    (Fraction(0, 1), Fraction(3, 1)),
                    (Fraction(3, 1), Fraction(0, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 4, 3, 2, 3, 2, 1, 2, 1, 0, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_triangle_case_2_2_4(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(6, 1), Fraction(6, 1)),
                    (Fraction(3, 1), Fraction(0, 1)),
                    (Fraction(0, 1), Fraction(3, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_triangle(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 4, 3, 2, 3, 2, 1, 2, 1, 0, 1]
        self.assertEqual(result, want)


class TestAxisArbitraryConvexCase(unittest.TestCase):
    def test_non_rectangular_convolution_convex_polygon_triangle(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(6, 1), Fraction(6, 1)),
                    (Fraction(3, 1), Fraction(0, 1)),
                    (Fraction(0, 1), Fraction(3, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_convex_polygon(
            list1, list2, geometry, prime)
        want = [0, 0, 0, 4, 3, 2, 3, 2, 1, 2, 1, 0, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_convex_polygon_quadrilateral(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(0, 1), Fraction(0, 1)),
                    (Fraction(4, 1), Fraction(2, 1)),
                    (Fraction(6, 1), Fraction(4, 1)),
                    (Fraction(2, 1), Fraction(4, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_convex_polygon(
            list1, list2, geometry, prime)
        want = [1, 0, 1, 2, 1, 2, 3, 2, 2, 1, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_convex_polygon_12_edges(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(3, 1), Fraction(0, 1)),
                    (Fraction(4, 1), Fraction(0, 1)),
                    (Fraction(6, 1), Fraction(1, 1)),
                    (Fraction(7, 1), Fraction(3, 1)),
                    (Fraction(7, 1), Fraction(4, 1)),
                    (Fraction(6, 1), Fraction(6, 1)),
                    (Fraction(4, 1), Fraction(7, 1)),
                    (Fraction(3, 1), Fraction(7, 1)),
                    (Fraction(1, 1), Fraction(6, 1)),
                    (Fraction(0, 1), Fraction(4, 1)),
                    (Fraction(0, 1), Fraction(3, 1)),
                    (Fraction(1, 1), Fraction(1, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_convex_polygon(
            list1, list2, geometry, prime)
        want = [0, 0, 1, 4, 5, 4, 5, 6, 5, 4, 5, 4, 1, 0, 0]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_convex_polygon_13_edges(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(Fraction(3, 1), Fraction(0, 1)),
                    (Fraction(4, 1), Fraction(0, 1)),
                    (Fraction(11, 2), Fraction(1, 1)),
                    (Fraction(6, 1), Fraction(3, 2)),
                    (Fraction(7, 1), Fraction(3, 1)),
                    (Fraction(7, 1), Fraction(4, 1)),
                    (Fraction(6, 1), Fraction(6, 1)),
                    (Fraction(4, 1), Fraction(7, 1)),
                    (Fraction(3, 1), Fraction(7, 1)),
                    (Fraction(1, 1), Fraction(6, 1)),
                    (Fraction(0, 1), Fraction(4, 1)),
                    (Fraction(0, 1), Fraction(3, 1)),
                    (Fraction(1, 1), Fraction(1, 1))]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_convex_polygon(
            list1, list2, geometry, prime)
        want = [0, 0, 1, 4, 5, 4, 5, 5, 5, 4, 5, 4, 1, 0, 0]
        self.assertEqual(result, want)


if __name__ == '__main__':
    unittest.main()
