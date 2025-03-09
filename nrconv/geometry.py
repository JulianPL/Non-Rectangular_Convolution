#!/usr/bin/python3
"""This module handles geometry calculation.
"""

from math import ceil, floor
from typing import List, Tuple
from fractions import Fraction

Point = Tuple[Fraction, Fraction]

def rectangle_inscribed_int(coordinates: List[Point]) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    """Retrieves bounds for the minimal and maximal integer coordinates of the given polygon"""
    (x_min, y_min), (x_max, y_max) = rectangle_inscribed(coordinates)
    return (ceil(x_min), ceil(y_min)), (floor(x_max), floor(y_max))

def rectangle_inscribed(coordinates: List[Point]) -> Tuple[Point, Point]:
    """Retrieves bounds for the minimal and maximal coordinates of the given polygon"""
    x_min = min(coordinate[0] for coordinate in coordinates)
    y_min = min(coordinate[1] for coordinate in coordinates)
    x_max = max(coordinate[0] for coordinate in coordinates)
    y_max = max(coordinate[1] for coordinate in coordinates)
    return (x_min, y_min), (x_max, y_max)

def closer_point(reference: Point, option_a: Point, option_b: Point) -> Point:
    """
    Returns the point closer to the reference.

    Args:
        reference (Point): The reference point (X).
        option_a (Point): First candidate point (Y).
        option_b (Point): Second candidate point (Z).

    Returns:
        Point: The closer point (option_a or option_b).
    """
    a_squared_distance = sum((ref - a) ** 2 for ref, a in zip(reference, option_a))
    b_squared_distance = sum((ref - b) ** 2 for ref, b in zip(reference, option_b))
    return option_a if a_squared_distance <= b_squared_distance else option_b

def opposing_rect_vertex(reference: Point, diag_start: Point, diag_end: Point) -> Point:
    """Determines the opposite rectangle vertex relative to the reference point.

    Given two diagonally opposite corners of an axis-aligned rectangle (`diag_start`, `diag_end`),
    this function determines which of the remaining two rectangle vertices is on the opposing site
    of the diagonal as the `reference` point.

    Args:
        reference (Point): The point to compare.
        diag_start (Point): One of the diagonally opposite rectangle corners.
        diag_end (Point): The other diagonally opposite rectangle corner.

    Returns:
        Point: The rectangle corner on the opposing side of the diagonal as `reference`.

    Example:
        >>> opposing_rect_vertex((Fraction(9, 1), Fraction(0, 1)), (Fraction(0, 1), Fraction(0, 1)), (Fraction(10, 1), Fraction(1, 1)))
        (Fraction(0, 1), Fraction(1, 1))
    """
    other_a, other_b = (diag_end[0], diag_start[1]), (diag_start[0], diag_end[1])
    weight_ax, weight_ay = abs(other_a[0] - diag_start[0]), abs(other_a[1] - diag_end[1])
    part_ax, part_ay = abs(other_a[0] - reference[0]), abs(other_a[1] - reference[1])
    return other_a if (part_ax / weight_ax) + (part_ay / weight_ay) > 1 else other_b