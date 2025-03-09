#!/usr/bin/python3

from fractions import Fraction

import unittest

import nrconv

class TestRectangleInscribed(unittest.TestCase):
    def test_rectangle_inscribed_int_single_point_int(self):
        geometry = [(Fraction(2, 1), Fraction(7, 1))]
        expected = ((2, 7), (2, 7))
        actual = nrconv.rectangle_inscribed_int(geometry)
        self.assertEqual(expected, actual)

    def test_rectangle_inscribed_int_single_point_rat(self):
        geometry = [(Fraction(5, 3), Fraction(7, 2))]
        expected = ((2, 4), (1, 3))
        actual = nrconv.rectangle_inscribed_int(geometry)
        self.assertEqual(expected, actual)

    def test_rectangle_inscribed_int_triangle_flat(self):
        geometry = [(Fraction(19, 20), Fraction(21, 20)), (Fraction(21, 20), Fraction(21, 20)),
                    (Fraction(119, 20), Fraction(179, 20))]
        expected = ((1, 2), (5, 8))
        actual = nrconv.rectangle_inscribed_int(geometry)
        self.assertEqual(expected, actual)

    def test_rectangle_inscribed_int_quadrilateral(self):
        geometry = [(Fraction(23, 10), Fraction(9, 10)), (Fraction(70, 10), Fraction(40, 10)),
                    (Fraction(45, 10), Fraction(69, 10)), (Fraction(11, 10), Fraction(36, 10))]
        expected = ((2, 1), (7, 6))
        actual = nrconv.rectangle_inscribed_int(geometry)
        self.assertEqual(expected, actual)

    def test_rectangle_inscribed_single_point_int(self):
        geometry = [(Fraction(2, 1), Fraction(7, 1))]
        expected = ((Fraction(2, 1), Fraction(7, 1)), (Fraction(2, 1), Fraction(7, 1)))
        actual = nrconv.rectangle_inscribed(geometry)
        self.assertEqual(expected, actual)

    def test_rectangle_inscribed_single_point_rat(self):
        geometry = [(Fraction(5, 3), Fraction(7, 2))]
        expected = ((Fraction(5, 3), Fraction(7, 2)), (Fraction(5, 3), Fraction(7, 2)))
        actual = nrconv.rectangle_inscribed(geometry)
        self.assertEqual(expected, actual)


    def test_rectangle_inscribed_triangle_flat(self):
        geometry = [(Fraction(19, 20), Fraction(21, 20)), (Fraction(21, 20), Fraction(21, 20)),
                    (Fraction(119, 20), Fraction(179, 20))]
        expected = ((Fraction(19, 20), Fraction(21, 20)), (Fraction(119, 20), Fraction(179, 20)))
        actual = nrconv.rectangle_inscribed(geometry)
        self.assertEqual(expected, actual)

    def test_rectangle_inscribed_quadrilateral(self):
        geometry = [(Fraction(23, 10), Fraction(9, 10)), (Fraction(70, 10), Fraction(40, 10)),
                    (Fraction(45, 10), Fraction(69, 10)), (Fraction(11, 10), Fraction(36, 10))]
        expected = ((Fraction(11, 10), Fraction(9, 10)), (Fraction(70, 10), Fraction(69, 10)))
        actual = nrconv.rectangle_inscribed(geometry)
        self.assertEqual(expected, actual)

class TestCloserPoint(unittest.TestCase):
    def test_closer_point(self):
        reference = (Fraction(30, 10), Fraction(30, 10))
        option_a = (Fraction(10, 10), Fraction(28, 10))
        option_b = (Fraction(29, 10), Fraction(11, 10))
        expected = (Fraction(29, 10), Fraction(11, 10))
        actual = nrconv.closer_point(reference, option_a, option_b)
        self.assertEqual(expected, actual)

class TestOpposingRectVertex(unittest.TestCase):
    def test_opposing_rect_vertex_simple(self):
        reference = (Fraction(50, 10), Fraction(24, 10))
        diag_start = (Fraction(10, 10), Fraction(20, 10))
        diag_end = (Fraction(65, 10), Fraction(39, 10))
        expected = (Fraction(10, 10), Fraction(39, 10))
        actual = nrconv.opposing_rect_vertex(reference, diag_start, diag_end)
        self.assertEqual(expected, actual)

    def test_opposing_rect_vertex_(self):
        reference = (Fraction(20, 10), Fraction(56, 10))
        diag_start = (Fraction(9, 10), Fraction(57, 10))
        diag_end = (Fraction(165, 10), Fraction(39, 10))
        expected = (Fraction(9, 10), Fraction(39, 10))
        actual = nrconv.opposing_rect_vertex(reference, diag_start, diag_end)
        self.assertEqual(expected, actual)
