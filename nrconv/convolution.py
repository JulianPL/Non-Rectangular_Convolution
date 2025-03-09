#!/usr/bin/python3
"""This module calculates convolutions with non-rectangular geometry.
"""

from dataclasses import dataclass
from math import ceil, floor
from typing import List, Tuple, Callable

import sympy

Point = Tuple[float, float]

@dataclass
class ConvolutionStep:
    geometry: List[Point]
    function: Callable[[List[int], List[int], List[Point], int], Tuple[List[int], int]]
    is_positive: bool  # True for addition, False for subtraction

def is_integer(number: float) -> bool:
    """Checks if number is an integer"""
    return number == ceil(number)

def rectangle_inscribed_int(coordinates: List[Point]) -> Tuple[int, int, int, int]:
    """Retrieves bounds for the minimal and maximal integer coordinates of the given polygon"""
    x_min, y_min, x_max, y_max = rectangle_inscribed(coordinates)
    return ceil(x_min), ceil(y_min), floor(x_max), floor(y_max)

def rectangle_inscribed(coordinates: List[Point]) -> Tuple[float, float, float, float]:
    """Retrieves bounds for the minimal and maximal coordinates of the given polygon"""
    x_min = min(coordinate[0] for coordinate in coordinates)
    y_min = min(coordinate[1] for coordinate in coordinates)
    x_max = max(coordinate[0] for coordinate in coordinates)
    y_max = max(coordinate[1] for coordinate in coordinates)
    return x_min, y_min, x_max, y_max

def retrieve_convolution_size(coordinates: List[Point]) -> Tuple[int, int]:
    """Retrieves bounds for the size and the starting index of the convolved sequence with given polygon"""
    x_min, y_min, x_max, y_max = rectangle_inscribed_int(coordinates)
    conv_min, conv_max = x_min + y_min, x_max + y_max
    if (x_min > x_max) or (y_min > y_max):
        conv_size = 0
    else:
        conv_size = conv_max - conv_min + 1
    return conv_size, conv_min

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
        >>> opposing_rect_vertex((9, 0), (0, 0), (10, 1))
        (0, 1)
    """
    other_a, other_b = (diag_end[0], diag_start[1]), (diag_start[0], diag_end[1])
    weight_ax, weight_ay = abs(other_a[0] - diag_start[0]), abs(other_a[1] - diag_end[1])
    part_ax, part_ay = abs(other_a[0] - reference[0]), abs(other_a[1] - reference[1])
    return other_a if (part_ax / weight_ax) + (part_ay / weight_ay) > 1 else other_b

def add_subslice(main: Tuple[List[int], int], sub: Tuple[List[int], int]) -> Tuple[List[int], int]:
    """Calculates the sum of two slices main and sub with offsets
    such that the sub-indices being in main"""

    out = main[0].copy()
    main_slice, main_start = main[0], main[1]
    sub_slice, sub_start = sub[0], sub[1]
    main_end, sub_end = main_start + len(main_slice), sub_start + len(sub_slice)

    if sub_start < main_start:
        raise IndexError(f"Can't add a slice starting with {sub_start} onto a slice {main_slice}")

    if sub_end > main_end:
        raise IndexError(f"Can't add a slice ending with {sub_end} onto a slice {main_end}")

    for index, value in enumerate(sub_slice):
        out[sub_start - main_start + index] += value
    return out, main_start

def sub_subslice(main: Tuple[List[int], int], sub: Tuple[List[int], int]) -> Tuple[List[int], int]:
    """Calculates the difference of two slices main and sub with offsets
    such that the sub-indices being in main"""

    inverse = [-elm for elm in sub[0]], sub[1]
    return add_subslice(main, inverse)

def add_convolution(
        list1: List[int], list2: List[int],
        conv: List[int], conv_min: int,
        steps: List[ConvolutionStep], ntt_prime: int
) -> Tuple[List[int], int]:
    """Applies a sequence of convolution steps.

    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        conv (List[int]): The current convolved sequence.
        conv_min (int): The offset of the convolved sequence.
        steps (List[ConvolutionStep]): The convolution steps.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.

        Second, the offset of the first index of the convolution.
    """

    for step in steps:
        conv_part, conv_part_min = step.function(list1, list2, step.geometry, ntt_prime)

        if step.is_positive:
            conv, _ = add_subslice((conv, conv_min), (conv_part, conv_part_min))
        else:
            conv, _ = sub_subslice((conv, conv_min), (conv_part, conv_part_min))

    return conv, conv_min

def non_rectangular_convolution_edge(list1: List[int], list2: List[int],
                                     geometry: List[Point],
                                     _ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 1: Single edges.

    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Point]): Two opposing vertices
            of the underlying edge.
        _ntt_prime (int): The prime for the number theoretic transform (unused since there is no "real" convolution).

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.

        Second, the offset of the first index of the convolution.
    """

    start, end = geometry[0], geometry[1]
    x_min, y_min, x_max, _ = rectangle_inscribed_int([start, end])
    conv_size, conv_min = retrieve_convolution_size([start, end])

    if conv_size == 0:
        return [], conv_min

    conv = [0] * conv_size

    # edge is vertical (since x_min and x_max are rounded, we have to use the actual x-coordinates)
    if start[0] == end[0]:
        # edge does not contain integer points
        if not is_integer(start[0]):
            return [], conv_min
        # x_min==x_max==start[0]
        # edge from y_min to y_max (of length conv_size)
        for index in range(conv_size):
            conv[index] = list1[x_min] * list2[y_min + index]
        return conv, conv_min

    # make start -> end increase in x-direction.
    if start[0] > end[0]:
        start, end = end, start

    # add all integer coordinates into the convolution slice
    x_diff = end[0] - start[0]
    y_diff = end[1] - start[1]
    for x_index in range(x_min, x_max + 1):
        y_index = start[1] + (x_index - start[0]) * y_diff / x_diff
        if is_integer(y_index):
            y_index = ceil(y_index)
            conv[x_index + y_index - conv_min] \
                = conv[x_index + y_index - conv_min] + list1[x_index] * list2[y_index]
    return conv, conv_min

def non_rectangular_convolution_rectangle(
        list1: List[int], list2: List[int], geometry: List[Tuple[float,
        float]],
        ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 2: Axis-aligned rectangles.
    All edges are included.

    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Point]): Two opposing vertices
            of the underlying rectangle.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.

        Second, the offset of the first index of the convolution.
    """

    start, end = geometry[0], geometry[1]
    x_min, y_min, x_max, y_max = rectangle_inscribed_int([start, end])
    conv_size, conv_min = retrieve_convolution_size([start, end])

    if conv_size == 0:
        return [], conv_min

    return sympy.discrete.convolutions.convolution_ntt(
        list1[x_min:x_max + 1], list2[y_min:y_max + 1],
        ntt_prime), conv_min

def non_rectangular_convolution_triangle_axis_aligned(
        list1: List[int], list2: List[int], geometry: List[Tuple[float,
        float]],
        ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 3: Axis-aligned triangles.
    All edges are included.

    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Point]): The three vertices defining the
            underlying triangle.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.

        Second, the offset of the first index of the convolution.
    """

    A, B, C = geometry[0], geometry[1], geometry[2]
    x_min, y_min, x_max, y_max = rectangle_inscribed([A, B, C])
    x_average, y_average = (x_min + x_max) / 2, (y_min + y_max) / 2
    conv_size, conv_min = retrieve_convolution_size([A, B, C])

    x_cathetus = A[0] + B[0] + C[0] - x_min - x_max
    y_cathetus = A[1] + B[1] + C[1] - y_min - y_max
    x_not_cathetus = x_min + x_max - x_cathetus
    y_not_cathetus = y_min + y_max - y_cathetus

    # Case 0: Degenerated triangle:
    if (x_min == x_max) or (y_min == y_max):
        return non_rectangular_convolution_edge(list1, list2, [(x_min, y_min), (x_max, y_max)], ntt_prime)

    # Case 1: Small triangle:
    if conv_size == 0:
        return [], conv_min
    if conv_size == 1:
        x_I, y_I = ceil(x_min), ceil(y_min)
        x_quotient = abs(x_I - x_cathetus) / (x_max - x_min)
        y_quotient = abs(y_I - y_cathetus) / (y_max - y_min)
        # the relative distances from the catheti have to be at most 1
        # for the point to be in the triangle
        if x_quotient + y_quotient > 1:
            return [], conv_min
        return [list1[x_I] * list2[y_I]], conv_min

    # Case 2: Large triangle:
    conv = [0] * conv_size

    # Define all steps in a single list
    steps = [
        ConvolutionStep(
            geometry=[(x_cathetus, y_cathetus), (x_average, y_average)],
            function=non_rectangular_convolution_rectangle,
            is_positive=True
        ),
        ConvolutionStep(
            geometry=[(x_average, y_average), (x_cathetus, y_average), (x_cathetus, y_not_cathetus)],
            function=non_rectangular_convolution_triangle_axis_aligned,
            is_positive=True
        ),
        ConvolutionStep(
            geometry=[(x_average, y_cathetus), (x_average, y_average)],
            function=non_rectangular_convolution_edge,
            is_positive=False
        ),
        ConvolutionStep(
            geometry=[(x_average, y_average), (x_average, y_cathetus), (x_not_cathetus, y_cathetus)],
            function=non_rectangular_convolution_triangle_axis_aligned,
            is_positive=True
        ),
        ConvolutionStep(
            geometry=[(x_cathetus, y_average), (x_average, y_average)],
            function=non_rectangular_convolution_edge,
            is_positive=False
        )
    ]

    # Apply all convolution steps in a structured way
    conv, conv_min = add_convolution(list1, list2, conv, conv_min, steps, ntt_prime)

    return conv, conv_min

def non_rectangular_convolution_triangle(
        list1: List[int], list2: List[int], geometry: List[Tuple[float,
        float]],
        ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 4: arbitrary Triangles.
    All edges are included.
    
    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Point]): The three vertices defining the
            underlying triangle.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.
        
        Second, the offset of the first index of the convolution.
    """
    A, B, C = geometry[0], geometry[1], geometry[2]
    x_min, y_min, x_max, y_max = rectangle_inscribed([A, B, C])
    conv_size, conv_min = retrieve_convolution_size([A, B, C])

    # Case 0.1: Degenerated triangle:
    if (x_min == x_max) or (y_min == y_max):
        return non_rectangular_convolution_edge(list1, list2, [(x_min, y_min), (x_max, y_max)], ntt_prime)

    conv = [0] * conv_size

    vertex_collisions = []
    vertex_non_collisions = []

    for vertex in [A, B, C]:
        if vertex[0] in {x_min, x_max} and vertex[1] in {y_min, y_max}:
            vertex_collisions.append(vertex)
        else:
            vertex_non_collisions.append(vertex)

    # Case 0.2: Axis-aligned triangle:
    if len(vertex_collisions) == 3:
        return non_rectangular_convolution_triangle_axis_aligned(list1, list2, geometry, ntt_prime)

    if len(vertex_collisions) == 2:
        # Lemma 12, Case 1: Two vertices are on opposing vertices of the surrounding rectangle.
        if (vertex_collisions[0][0] != vertex_collisions[1][0]) and (
                vertex_collisions[0][1] != vertex_collisions[1][1]):
            # Compute the key geometry components dynamically
            quart = opposing_rect_vertex(vertex_non_collisions[0], vertex_collisions[0], vertex_collisions[1])
            base_points = [(vertex_non_collisions[0][0], quart[1]), (quart[0], vertex_non_collisions[0][1])]

            # Ensure base_point[i] is closest to vertex_collisions[i]
            base_points = [closer_point(vc, base_points[0], base_points[1]) for vc in vertex_collisions]

            # Define convolution operations
            add_rectangle = ConvolutionStep(
                geometry=[vertex_non_collisions[0], quart],
                function=non_rectangular_convolution_rectangle,
                is_positive=True)
            add_first_triangle = ConvolutionStep(
                geometry=[vertex_collisions[0], base_points[0], vertex_non_collisions[0]],
                function=non_rectangular_convolution_triangle_axis_aligned,
                is_positive=True)
            sub_first_edge = ConvolutionStep(
                geometry=[base_points[0], vertex_non_collisions[0]],
                function=non_rectangular_convolution_edge,
                is_positive=False)
            add_second_triangle = ConvolutionStep(
                geometry=[vertex_collisions[1], base_points[1], vertex_non_collisions[0]],
                function=non_rectangular_convolution_triangle_axis_aligned,
                is_positive=True)
            sub_second_edge = ConvolutionStep(
                geometry=[base_points[1], vertex_non_collisions[0]],
                function=non_rectangular_convolution_edge,
                is_positive=False)
            sub_third_triangle = ConvolutionStep(
                geometry=[vertex_collisions[0], vertex_collisions[1], quart],
                function=non_rectangular_convolution_triangle_axis_aligned,
                is_positive=False)
            add_hypotenuse_edge = ConvolutionStep(
                geometry=[vertex_collisions[0], vertex_collisions[1]],
                function=non_rectangular_convolution_edge,
                is_positive=True)

            # Store in steps list
            steps = [add_rectangle,
                     add_first_triangle, sub_first_edge,
                     add_second_triangle, sub_second_edge,
                     sub_third_triangle, add_hypotenuse_edge]

        # Case 2.1: Two vertices are on the same edge of the surrounding rectangle.
        else:
            # calculate base point on the common edge
            if vertex_collisions[0][0] == vertex_collisions[1][0]:
                base_point = (vertex_collisions[0][0], vertex_non_collisions[0][1])
            else:
                base_point = (vertex_non_collisions[0][0], vertex_collisions[0][1])

            add_first_triangle = ConvolutionStep(
                geometry=[vertex_collisions[0], base_point, vertex_non_collisions[0]],
                function=non_rectangular_convolution_triangle_axis_aligned,
                is_positive=True)
            add_second_triangle = ConvolutionStep(
                geometry=[vertex_collisions[1], base_point, vertex_non_collisions[0]],
                function=non_rectangular_convolution_triangle_axis_aligned,
                is_positive=True)
            sub_common_edge = ConvolutionStep(
                geometry=[base_point, vertex_non_collisions[0]],
                function=non_rectangular_convolution_edge,
                is_positive=False)

            steps = [add_first_triangle, add_second_triangle, sub_common_edge]

        conv, conv_min = add_convolution(list1, list2, conv, conv_min, steps, ntt_prime)

        return conv, conv_min

    # Case 2.2: Only one vertex of the triangle coincides with the surrounding rectangle.
    opposite_collision = (x_min + x_max - vertex_collisions[0][0], y_min + y_max - vertex_collisions[0][1])

    corner0 = opposing_rect_vertex(vertex_non_collisions[1], vertex_collisions[0], opposite_collision)
    corner1 = opposing_rect_vertex(vertex_non_collisions[0], vertex_collisions[0], opposite_collision)

    add_rectangle = ConvolutionStep(
        geometry=[vertex_collisions[0], opposite_collision],
        function=non_rectangular_convolution_rectangle,
        is_positive=True)
    sub_first_triangle = ConvolutionStep(
        geometry=[vertex_collisions[0], corner0, vertex_non_collisions[0]],
        function=non_rectangular_convolution_triangle_axis_aligned,
        is_positive=False)
    add_first_edge = ConvolutionStep(
        geometry=[vertex_collisions[0], vertex_non_collisions[0]],
        function=non_rectangular_convolution_edge,
        is_positive=True)
    sub_second_triangle = ConvolutionStep(
        geometry=[vertex_collisions[0], corner1, vertex_non_collisions[1]],
        function=non_rectangular_convolution_triangle_axis_aligned,
        is_positive=False)
    add_second_edge = ConvolutionStep(
        geometry=[vertex_collisions[0], vertex_non_collisions[1]],
        function=non_rectangular_convolution_edge,
        is_positive=True)
    sub_third_triangle = ConvolutionStep(
        geometry=[opposite_collision, vertex_non_collisions[0], vertex_non_collisions[1]],
        function=non_rectangular_convolution_triangle_axis_aligned,
        is_positive=False)
    add_third_edge = ConvolutionStep(
        geometry=[vertex_non_collisions[0], vertex_non_collisions[1]],
        function=non_rectangular_convolution_edge,
        is_positive=True)

    steps = [add_rectangle,
             sub_first_triangle, add_first_edge,
             sub_second_triangle, add_second_edge,
             sub_third_triangle, add_third_edge]

    conv, conv_min = add_convolution(list1, list2, conv, conv_min, steps, ntt_prime)

    return conv, conv_min

def non_rectangular_convolution_convex_polygon(
        list1: List[int], list2: List[int], geometry: List[Tuple[float,
        float]],
        ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 5: Arbitrary convex polygons.
    All edges are included.
    
    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Point]): The vertices defining the
            underlying polygon.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.
        
        Second, the offset of the first index of the convolution.
    """
    number_vertices = len(geometry)

    if number_vertices == 2:
        return non_rectangular_convolution_edge(list1, list2, geometry, ntt_prime)

    if number_vertices == 3:
        return non_rectangular_convolution_triangle(list1, list2, geometry, ntt_prime)

    conv_size, conv_min = retrieve_convolution_size(geometry)

    conv = [0] * conv_size

    if number_vertices == 4:
        # Add first triangle.
        # Add second triangle.
        # Subtract edge between the two triangles, which was counted twice.
        steps = [
            ConvolutionStep(
                geometry=[geometry[0], geometry[1], geometry[2]],
                function=non_rectangular_convolution_triangle,
                is_positive=True
            ),
            ConvolutionStep(
                geometry=[geometry[2], geometry[3], geometry[0]],
                function=non_rectangular_convolution_triangle,
                is_positive=True
            ),
            ConvolutionStep(
                geometry=[geometry[0], geometry[2]],
                function=non_rectangular_convolution_edge,
                is_positive=False
            )
        ]

        conv, conv_min = add_convolution(list1, list2, conv, conv_min, steps, ntt_prime)

        return conv, conv_min

    # Polygon v_{0} - v{2} - v{4} - v_{6} - ... - v_{2*floor(k/2)}:
    steps = [
        ConvolutionStep(
            geometry=geometry[::2],
            function=non_rectangular_convolution_convex_polygon,
            is_positive=True
        )
    ]
    for index in range(0, number_vertices - 2, 2):
        steps.append(ConvolutionStep(
            geometry=[geometry[index], geometry[index + 1], geometry[index + 2]],
            function=non_rectangular_convolution_triangle,
            is_positive=True
        ))
        steps.append(ConvolutionStep(
            geometry=[geometry[index], geometry[index + 2]],
            function=non_rectangular_convolution_edge,
            is_positive=False
        ))
    if number_vertices % 2 == 0:
        steps.append(ConvolutionStep(
            geometry=[geometry[number_vertices - 2], geometry[number_vertices - 1], geometry[0]],
            function=non_rectangular_convolution_triangle,
            is_positive=True
        ))
        steps.append(ConvolutionStep(
            geometry=[geometry[number_vertices - 2], geometry[0]],
            function=non_rectangular_convolution_edge,
            is_positive=False
        ))

    conv, conv_min = add_convolution(list1, list2, conv, conv_min, steps, ntt_prime)

    return conv, conv_min
