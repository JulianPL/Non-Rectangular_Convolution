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
    """Non-Rectangular Convolution -- Base Case 1: single Edge.
    
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
    """Non-Rectangular Convolution -- Base Case 2: axis-aligned Rectangle.
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
    """Non-Rectangular Convolution -- Base Case 3: axis-aligned Triangle.
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
        return non_rectangular_convolution_rectangle(list1, list2,
                                                     ((x_min, y_min),
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

    # add rectangle
    conv_rectangle, conv_rectangle_min = non_rectangular_convolution_rectangle(
        list1, list2, ((x_cathetus, y_cathetus), (x_average, y_average)),
        ntt_prime)
    for index, value in enumerate(conv_rectangle):
        conv[conv_rectangle_min - conv_min +
             index] = conv[conv_rectangle_min - conv_min + index] + value

    # add first triangle
    conv_triangle1, conv_triangle1_min = non_rectangular_convolution_triangle_axis_aligned(
        list1, list2, ((x_average, y_average), (x_cathetus, y_average),
                       (x_cathetus, y_not_cathetus)), ntt_prime)
    for index, value in enumerate(conv_triangle1):
        conv[conv_triangle1_min - conv_min +
             index] = conv[conv_triangle1_min - conv_min + index] + value

    # subtract line between rectangle and first triangle, which was counted twice.
    conv_line1, conv_line1_min = non_rectangular_convolution_rectangle(
        list1, list2, ((x_average, y_cathetus), (x_average, y_average)),
        ntt_prime)
    for index, value in enumerate(conv_line1):
        conv[conv_line1_min - conv_min +
             index] = conv[conv_line1_min - conv_min + index] - value

    # add second triangle
    conv_triangle2, conv_triangle2_min = non_rectangular_convolution_triangle_axis_aligned(
        list1, list2, ((x_average, y_average), (x_average, y_cathetus),
                       (x_not_cathetus, y_cathetus)), ntt_prime)
    for index, value in enumerate(conv_triangle2):
        conv[conv_triangle2_min - conv_min +
             index] = conv[conv_triangle2_min - conv_min + index] + value

    # subtract line between rectangle and second triangle, which was counted twice.
    conv_line2, conv_line2_min = non_rectangular_convolution_rectangle(
        list1, list2, ((x_cathetus, y_average), (x_average, y_average)),
        ntt_prime)
    for index, value in enumerate(conv_line2):
        conv[conv_line2_min - conv_min +
             index] = conv[conv_line2_min - conv_min + index] - value

    return conv, conv_min
