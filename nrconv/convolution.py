#!/usr/bin/python3
"""This module calculates convolutions with non-rectangular geometry.
"""

from typing import List, Tuple

import sympy

from math import ceil, floor


def is_integer(number: float) -> bool:
    """Checks if number is an integer"""
    return number == ceil(number)


def non_rectangular_convolution_edge(list1: List[int], list2: List[int],
                                     geometry: List[Tuple[float, float]],
                                     ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 1: single Edges.
    
    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Tuple[float, float]]): Two opposing vertices 
            of the underlying edge.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.
        
        Second, the offset of the first index of the convolution.
    """

    A = geometry[0]
    B = geometry[1]

    x_min = ceil(min(A[0], B[0]))
    y_min = ceil(min(A[1], B[1]))
    x_max = floor(max(A[0], B[0]))
    y_max = floor(max(A[1], B[1]))

    if (x_min > x_max) or (y_min > y_max):
        return [], x_min + y_min

    conv_min = x_min + y_min
    conv_max = x_max + y_max
    conv_size = conv_max - conv_min + 1
    conv = [0] * conv_size

    if A[0] == B[0]:
        if is_integer(A[0]):
            for index in range(conv_size):
                conv[index] = list1[ceil(A[0])] * list2[y_min + index]
            return conv, conv_min
        else:
            return [], conv_min

    # make A -> B increase in x-direction.
    if A[0] > B[0]:
        temp_vertex = A
        A = B
        B = temp_vertex

    x_diff = B[0] - A[0]
    y_diff = B[1] - A[1]
    for x_index in range(x_min, x_max + 1):
        y_index = A[1] + (x_index - A[0]) * y_diff / x_diff
        if is_integer(y_index):
            y_index = ceil(y_index)
            conv[x_index + y_index -
                 conv_min] = conv[x_index + y_index -
                                  conv_min] + list1[x_index] * list2[y_index]

    return conv, conv_min


def non_rectangular_convolution_rectangle(
        list1: List[int], list2: List[int], geometry: List[Tuple[float,
                                                                 float]],
        ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 2: axis-aligned Rectangles.
    All edges are included.
    
    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Tuple[float, float]]): Two opposing vertices 
            of the underlying rectangle.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.
        
        Second, the offset of the first index of the convolution.
    """

    A = geometry[0]
    B = geometry[1]

    x_min = ceil(min(A[0], B[0]))
    y_min = ceil(min(A[1], B[1]))
    x_max = floor(max(A[0], B[0]))
    y_max = floor(max(A[1], B[1]))

    if (x_min > x_max) or (y_min > y_max):
        return [], x_min + y_min

    return sympy.discrete.convolutions.convolution_ntt(
        list1[x_min:x_max + 1], list2[y_min:y_max + 1],
        ntt_prime), x_min + y_min


def non_rectangular_convolution_triangle_axis_aligned(
        list1: List[int], list2: List[int], geometry: List[Tuple[float,
                                                                 float]],
        ntt_prime: int) -> Tuple[List[int], int]:
    """Non-Rectangular Convolution -- Base Case 3: axis-aligned Triangles.
    All edges are included.
    
    Args:
        list1 (List[int]): The first list.
        list2 (List[int]): The second list.
        geometry (List[Tuple[float, float]]): The three vertices defining the
            underlying triangle.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.
        
        Second, the offset of the first index of the convolution.
    """

    A = geometry[0]
    B = geometry[1]
    C = geometry[2]

    x_min = min(A[0], B[0], C[0])
    y_min = min(A[1], B[1], C[1])
    x_max = max(A[0], B[0], C[0])
    y_max = max(A[1], B[1], C[1])
    x_average = (x_min + x_max) / 2
    y_average = (y_min + y_max) / 2
    x_cathetus = A[0] + B[0] + C[0] - x_min - x_max
    y_cathetus = A[1] + B[1] + C[1] - y_min - y_max
    x_not_cathetus = x_min + x_max - x_cathetus
    y_not_cathetus = y_min + y_max - y_cathetus

    # Case 0: Degenerated triangle:
    if (x_min == x_max) or (y_min == y_max):
        return non_rectangular_convolution_edge(list1, list2, ((x_min, y_min),
                                                               (x_max, y_max)),
                                                ntt_prime)

    # Case 1: Small triangle:
    if (x_max - x_min) + (y_max - y_min) < 1:
        x_I = ceil(x_min)
        y_I = ceil(y_min)

        if (x_I > x_max) or (y_I > y_max):
            return [], x_I + y_I

        if (abs(x_I - x_cathetus) / (x_max - x_min)) + (abs(y_I - y_cathetus) /
                                                        (y_max - y_min)) > 1:
            return [], x_I + y_I

        return [list1[x_I] * list2[y_I]], x_I + y_I

    # Case 2: Large triangle:
    conv_min = ceil(x_min) + ceil(y_min)
    conv_max = floor(x_max) + floor(y_max)
    conv_size = conv_max - conv_min + 1
    conv = [0] * conv_size

    # add rectangle.
    conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
        list1, list2, ((x_cathetus, y_cathetus), (x_average, y_average)),
        ntt_prime)
    for index, value in enumerate(conv_rectangle):
        conv[conv_rectangle_min - conv_min +
             index] = conv[conv_rectangle_min - conv_min + index] + value

    # add first triangle.
    conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
        list1, list2, ((x_average, y_average), (x_cathetus, y_average),
                       (x_cathetus, y_not_cathetus)), ntt_prime)
    for index, value in enumerate(conv_triangle1):
        conv[conv_triangle1_min - conv_min +
             index] = conv[conv_triangle1_min - conv_min + index] + value

    # subtract edge between rectangle and first triangle, which was counted twice.
    conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
        list1, list2, ((x_average, y_cathetus), (x_average, y_average)),
        ntt_prime)
    for index, value in enumerate(conv_edge1):
        conv[conv_edge1_min - conv_min +
             index] = conv[conv_edge1_min - conv_min + index] - value

    # add second triangle.
    conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
        list1, list2, ((x_average, y_average), (x_average, y_cathetus),
                       (x_not_cathetus, y_cathetus)), ntt_prime)
    for index, value in enumerate(conv_triangle2):
        conv[conv_triangle2_min - conv_min +
             index] = conv[conv_triangle2_min - conv_min + index] + value

    # subtract edge between rectangle and second triangle, which was counted twice.
    conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
        list1, list2, ((x_cathetus, y_average), (x_average, y_average)),
        ntt_prime)
    for index, value in enumerate(conv_edge2):
        conv[conv_edge2_min - conv_min +
             index] = conv[conv_edge2_min - conv_min + index] - value

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
        geometry (List[Tuple[float, float]]): The three vertices defining the
            underlying triangle.
        ntt_prime (int): The prime for the number theoretic transform.

    Returns:
        First, the convolution of the two lists with the given
            base geometry as a list of integers.
        
        Second, the offset of the first index of the convolution.
    """

    A = geometry[0]
    B = geometry[1]
    C = geometry[2]

    x_min = min(A[0], B[0], C[0])
    y_min = min(A[1], B[1], C[1])
    x_max = max(A[0], B[0], C[0])
    y_max = max(A[1], B[1], C[1])

    # Case 0.1: degenerated triangle:
    if (x_min == x_max) or (y_min == y_max):
        return non_rectangular_convolution_edge(list1, list2, ((x_min, y_min),
                                                               (x_max, y_max)),
                                                ntt_prime)

    conv_min = ceil(x_min) + ceil(y_min)
    conv_max = floor(x_max) + floor(y_max)
    conv_size = conv_max - conv_min + 1
    conv = [0] * conv_size

    vertex_collisions = []
    vertex_non_collisions = []
    for vertex in [A, B, C]:
        collision = False
        for x_extreme in [x_min, x_max]:
            for y_extreme in [y_min, y_max]:
                if vertex == (x_extreme, y_extreme):
                    vertex_collisions.append(vertex)
                    collision = True
        if not collision:
            vertex_non_collisions.append(vertex)

    # Case 0.2: axis-aligned triangle:
    if len(vertex_collisions) == 3:
        return non_rectangular_convolution_triangle_axis_aligned(
            list1, list2, geometry, ntt_prime)

    if len(vertex_collisions) == 2:
        # Lemma 12, Case 1: two vertices are on opposing vertices of the surrounding rectangle.
        if (vertex_collisions[0][0] != vertex_collisions[1][0]) and (
                vertex_collisions[0][1] != vertex_collisions[1][1]):
            # Case 1.1: x_min/y_min -- x_max/y_max.
            if (x_min == vertex_collisions[0][0]
                    and y_min == vertex_collisions[0][1]) or (
                        x_max == vertex_collisions[0][0]
                        and y_max == vertex_collisions[0][1]):
                # Case 1.1.1: third vertex in the bottom right.
                if (abs(vertex_non_collisions[0][0] - x_min) /
                    (x_max - x_min)) + (
                        abs(vertex_non_collisions[0][1] - y_max) /
                        (y_max - y_min)) > 1:
                    conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
                        list1, list2,
                        (vertex_non_collisions[0], (x_min, y_max)), ntt_prime)
                    conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_min, vertex_non_collisions[0][1]),
                         (x_min, y_min), vertex_non_collisions[0]), ntt_prime)
                    conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_min, vertex_non_collisions[0][1]),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_max, y_max), (vertex_non_collisions[0][0], y_max),
                         vertex_non_collisions[0]), ntt_prime)
                    conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
                        list1, list2, ((vertex_non_collisions[0][0], y_max),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle3, conv_triangle3_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_min, y_max), (x_min, y_min), (x_max, y_max)),
                        ntt_prime)
                    conv_edge3, conv_edge3_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_min, y_min), (x_max, y_max)),
                        ntt_prime)
                # Case 1.1.2: third vertex in the top left.
                else:
                    conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
                        list1, list2,
                        (vertex_non_collisions[0], (x_max, y_min)), ntt_prime)
                    conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_max, vertex_non_collisions[0][1]),
                         (x_max, y_max), vertex_non_collisions[0]), ntt_prime)
                    conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_max, vertex_non_collisions[0][1]),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_min, y_min), (vertex_non_collisions[0][0], y_min),
                         vertex_non_collisions[0]), ntt_prime)
                    conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
                        list1, list2, ((vertex_non_collisions[0][0], y_min),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle3, conv_triangle3_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_max, y_min), (x_min, y_min), (x_max, y_max)),
                        ntt_prime)
                    conv_edge3, conv_edge3_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_min, y_min), (x_max, y_max)),
                        ntt_prime)
            # Case 1.2: x_min/y_min -- x_max/y_max.
            if (x_min == vertex_collisions[0][0]
                    and y_max == vertex_collisions[0][1]) or (
                        x_max == vertex_collisions[0][0]
                        and y_min == vertex_collisions[0][1]):
                # Case 1.2.1: third vertex in the top right.
                if (abs(vertex_non_collisions[0][0] - x_min) /
                    (x_max - x_min)) + (
                        abs(vertex_non_collisions[0][1] - y_min) /
                        (y_max - y_min)) > 1:
                    conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
                        list1, list2,
                        (vertex_non_collisions[0], (x_min, y_min)), ntt_prime)
                    conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_min, vertex_non_collisions[0][1]),
                         (x_min, y_max), vertex_non_collisions[0]), ntt_prime)
                    conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_min, vertex_non_collisions[0][1]),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_max, y_min), (vertex_non_collisions[0][0], y_min),
                         vertex_non_collisions[0]), ntt_prime)
                    conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
                        list1, list2, ((vertex_non_collisions[0][0], y_min),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle3, conv_triangle3_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_min, y_max), (x_min, y_min), (x_max, y_min)),
                        ntt_prime)
                    conv_edge3, conv_edge3_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_min, y_max), (x_max, y_min)),
                        ntt_prime)
                # Case 1.2.2: third vertex in the bottom left.
                else:
                    conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
                        list1, list2,
                        (vertex_non_collisions[0], (x_max, y_max)), ntt_prime)
                    conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_max, vertex_non_collisions[0][1]),
                         (x_max, y_min), vertex_non_collisions[0]), ntt_prime)
                    conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_max, vertex_non_collisions[0][1]),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_min, y_max), (vertex_non_collisions[0][0], y_max),
                         vertex_non_collisions[0]), ntt_prime)
                    conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
                        list1, list2, ((vertex_non_collisions[0][0], y_max),
                                       vertex_non_collisions[0]), ntt_prime)
                    conv_triangle3, conv_triangle3_min = non_rectangular_convolution_triangle_axis_aligned(
                        list1, list2,
                        ((x_min, y_max), (x_max, y_max), (x_max, y_min)),
                        ntt_prime)
                    conv_edge3, conv_edge3_min = non_rectangular_convolution_edge(
                        list1, list2, ((x_min, y_max), (x_max, y_min)),
                        ntt_prime)

            # Summing the convolutions for Case 1:
            # add rectangle.
            for index, value in enumerate(conv_rectangle):
                conv[conv_rectangle_min - conv_min +
                     index] = conv[conv_rectangle_min - conv_min +
                                   index] + value
            # add first triangle.
            for index, value in enumerate(conv_triangle1):
                conv[conv_triangle1_min - conv_min +
                     index] = conv[conv_triangle1_min - conv_min +
                                   index] + value
            # subtract edge between rectangle and first triangle, which was counted twice.
            for index, value in enumerate(conv_edge1):
                conv[conv_edge1_min - conv_min +
                     index] = conv[conv_edge1_min - conv_min + index] - value
            # add second triangle.
            for index, value in enumerate(conv_triangle2):
                conv[conv_triangle2_min - conv_min +
                     index] = conv[conv_triangle2_min - conv_min +
                                   index] + value
            # subtract edge between rectangle and second triangle, which was counted twice.
            for index, value in enumerate(conv_edge2):
                conv[conv_edge2_min - conv_min +
                     index] = conv[conv_edge2_min - conv_min + index] - value
            # subtract third triangle.
            for index, value in enumerate(conv_triangle3):
                conv[conv_triangle3_min - conv_min +
                     index] = conv[conv_triangle3_min - conv_min +
                                   index] - value
            # add hypothenuse of the third triangle, which was wrongly subtracted.
            for index, value in enumerate(conv_edge3):
                conv[conv_edge3_min - conv_min +
                     index] = conv[conv_edge3_min - conv_min + index] + value
            return conv, conv_min

        # Case 2.1: two vertices are on the same edge of the surrounding rectangle.
        else:
            conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
                list1, list2, ((x_min, y_min), (x_max, y_max)), ntt_prime)
            # Case 2.1.1 vertical edge.
            if vertex_collisions[0][0] == vertex_collisions[1][0]:
                conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
                    list1, list2,
                    ((x_min, y_min), (x_max, y_min), vertex_non_collisions[0]),
                    ntt_prime)
                conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
                    list1, list2,
                    ((x_min + x_max - vertex_non_collisions[0][0], y_min),
                     vertex_non_collisions[0]), ntt_prime)
                conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
                    list1, list2,
                    ((x_min, y_max), (x_max, y_max), vertex_non_collisions[0]),
                    ntt_prime)
                conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
                    list1, list2,
                    ((x_min + x_max - vertex_non_collisions[0][0], y_max),
                     vertex_non_collisions[0]), ntt_prime)
            # Case 2.1.2 horizontal edge.
            else:
                conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
                    list1, list2,
                    ((x_min, y_min), (x_min, y_max), vertex_non_collisions[0]),
                    ntt_prime)
                conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
                    list1, list2,
                    ((x_min, y_min + y_max - vertex_non_collisions[0][1]),
                     vertex_non_collisions[0]), ntt_prime)
                conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
                    list1, list2,
                    ((x_max, y_min), (x_max, y_max), vertex_non_collisions[0]),
                    ntt_prime)
                conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
                    list1, list2,
                    ((x_max, y_min + y_max - vertex_non_collisions[0][1]),
                     vertex_non_collisions[0]), ntt_prime)

            # Summing the convolutions for Case 2.1:
            # add rectangle.
            for index, value in enumerate(conv_rectangle):
                conv[conv_rectangle_min - conv_min +
                     index] = conv[conv_rectangle_min - conv_min +
                                   index] + value
            # subtract first triangle.
            for index, value in enumerate(conv_triangle1):
                conv[conv_triangle1_min - conv_min +
                     index] = conv[conv_triangle1_min - conv_min +
                                   index] - value
            # add edge between given triange and first triangle, which was wrongly subtracted.
            for index, value in enumerate(conv_edge1):
                conv[conv_edge1_min - conv_min +
                     index] = conv[conv_edge1_min - conv_min + index] + value
            # subtract second triangle.
            for index, value in enumerate(conv_triangle2):
                conv[conv_triangle2_min - conv_min +
                     index] = conv[conv_triangle2_min - conv_min +
                                   index] - value
            # add edge between given triange and second triangle, which was wrongly subtracted.
            for index, value in enumerate(conv_edge2):
                conv[conv_edge2_min - conv_min +
                     index] = conv[conv_edge2_min - conv_min + index] + value
            return conv, conv_min

    # Case 2.2: only one vertex of the triangle coincides with the surrounding rectangle.
    # Switch non_collision-vertices such that the first is horizontally opposed to verex_collisions.
    if (vertex_non_collisions[0][0] !=
            x_max + x_min - vertex_collisions[0][0]):
        temp = vertex_non_collisions[0]
        vertex_non_collisions[0] = vertex_non_collisions[1]
        vertex_non_collisions[1] = temp

    conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
        list1, list2, ((x_min, y_min), (x_max, y_max)), ntt_prime)
    conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
        list1, list2, (vertex_collisions[0],
                       (vertex_non_collisions[0][0], vertex_collisions[0][1]),
                       vertex_non_collisions[0]), ntt_prime)
    conv_edge1, conv_edge1_min = non_rectangular_convolution_edge(
        list1, list2, (vertex_collisions[0], vertex_non_collisions[0]),
        ntt_prime)
    conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
        list1, list2, (vertex_collisions[0],
                       (vertex_collisions[0][0], vertex_non_collisions[1][1]),
                       vertex_non_collisions[1]), ntt_prime)
    conv_edge2, conv_edge2_min = non_rectangular_convolution_edge(
        list1, list2, (vertex_collisions[0], vertex_non_collisions[1]),
        ntt_prime)
    conv_triangle3, conv_triangle3_min = non_rectangular_convolution_triangle_axis_aligned(
        list1, list2, ((x_min + x_max - vertex_collisions[0][0],
                        y_min + y_max - vertex_collisions[0][1]),
                       vertex_non_collisions[0], vertex_non_collisions[1]),
        ntt_prime)
    conv_edge3, conv_edge3_min = non_rectangular_convolution_edge(
        list1, list2, (vertex_non_collisions[0], vertex_non_collisions[1]),
        ntt_prime)

    # Summing the convolutions for Case 2.2:
    # add rectangle.
    for index, value in enumerate(conv_rectangle):
        conv[conv_rectangle_min - conv_min +
             index] = conv[conv_rectangle_min - conv_min + index] + value
    # subtract first triangle.
    for index, value in enumerate(conv_triangle1):
        conv[conv_triangle1_min - conv_min +
             index] = conv[conv_triangle1_min - conv_min + index] - value
    # add edge between given triange and first triangle, which was wrongly subtracted.
    for index, value in enumerate(conv_edge1):
        conv[conv_edge1_min - conv_min +
             index] = conv[conv_edge1_min - conv_min + index] + value
    # subtract second triangle.
    for index, value in enumerate(conv_triangle2):
        conv[conv_triangle2_min - conv_min +
             index] = conv[conv_triangle2_min - conv_min + index] - value
    # add edge between given triange and second triangle, which was wrongly subtracted.
    for index, value in enumerate(conv_edge2):
        conv[conv_edge2_min - conv_min +
             index] = conv[conv_edge2_min - conv_min + index] + value
    # subtract third triangle.
    for index, value in enumerate(conv_triangle3):
        conv[conv_triangle3_min - conv_min +
             index] = conv[conv_triangle3_min - conv_min + index] - value
    # add edge between given triange and third triangle, which was wrongly subtracted.
    for index, value in enumerate(conv_edge3):
        conv[conv_edge3_min - conv_min +
             index] = conv[conv_edge3_min - conv_min + index] + value
    return conv, conv_min
