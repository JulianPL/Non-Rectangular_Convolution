#!/usr/bin/python3
"""This module calculates convolutions with non-rectangular geometry.
"""

from typing import List, Tuple

import sympy

from math import ceil, floor


def non_rectangular_convolution_rectangle(list1: List[int], list2: List[int],
                                          geometry: List[Tuple[float, float]],
                                          ntt_prime: int):
    """Non-Rectangular Convolution -- Base Case: axis-aligned Rectangle.
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
