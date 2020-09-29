#!/usr/bin/python3

from fractions import Fraction

import unittest

import nrconv

import sympy


class TestRectangleCase(unittest.TestCase):
    def test_non_rectangular_convolution_rectangle_full(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(0, 0), (7, 7)]
        prime = nrconv.create_ntt_prime(list1, list2)
        result, _ = nrconv.convolution.non_rectangular_convolution_rectangle(
            list1, list2, geometry, prime)
        want = [1, 2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 1]
        self.assertEqual(result, want)

    def test_non_rectangular_convolution_rectangle_int_part(self):
        list1 = [1, 1, 1, 1, 1, 1, 1, 1]
        list2 = [1, 1, 1, 1, 1, 1, 1, 1]
        geometry = [(1, 2), (4, 6)]
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
        geometry = [(2, 5), (2, 5)]
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


if __name__ == '__main__':
    unittest.main()
