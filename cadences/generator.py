#!/usr/bin/python3
"""This module calculates lists whose underlying string have or avoid certain cadence patterns.
"""

from random import randrange
from typing import List

import matplotlib.pyplot as plt


def cadence_index_plausibility(start: int, distance: int, length_cadence: int, length_list: int) -> bool:
    # index plausibility check
    if start - distance >= 0:
        return False
    if start < 0:
        return False
    if start + (length_cadence - 1) * distance >= length_list:
        return False
    if start + length_cadence * distance < length_list:
        return False
    return True


def check_cadence(list_string: List[int], char: int, start: int, distance: int, length_cadence: int) -> bool:
    if not cadence_index_plausibility(start, distance, length_cadence, len(list_string)):
        return False
    # check characters
    for i in range(start, start + (length_cadence - 1) * distance, distance):
        if list_string[i] != char:
            return False
    return True


def check_partial_cadence(list_string: List[int], start: int, distance: int, length_cadence: int, partial: List[int]) -> bool:
    if not cadence_index_plausibility(start, distance, length_cadence, len(list_string)):
        return False
    # check characters
    for i in partial:
        if list_string[start + i * distance] == 0:
            return False
    return True


def make_4_cadence_free(min_length: int) -> List[int]:
    length_list = min_length if min_length % 7 == 0 else (1 + (min_length // 7)) * 7
    length_block = length_list // 7
    ret = [1 if (i // length_block) % 2 == 0 else 0 for i in range(length_list)]
    for distance_raw in range(length_block):
        for start in range(length_block):
            for distance in [2 * length_block + distance_raw, 2 * length_block - distance_raw]:
                if check_cadence(ret, 1, start, distance, 4):
                    ret[start + randrange(4) * distance] = 0
    return ret


def reduce_non_partial_4_cadences(base: List[int]) -> List[int]:
    length_block = len(base) // 7
    for distance_raw in range(length_block):
        for start in range(length_block):
            for distance in [2 * length_block + distance_raw, 2 * length_block - distance_raw]:
                for partial_num in range(4):
                    partial = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]][partial_num]
                    if check_partial_cadence(base, start, distance, 4, partial):
                        for i in partial:
                            base[start + i * distance] = base[start + i * distance] | (1 << (partial_num + 1))
    for i in range(len(base)):
        if bin(base[i]).count('1') == 4:
            base[i] = 1
        else:
            base[i] = 0
    return base

def test():
    a = []
    b = []
    c = []
    d = []
    e = []
    for i in range(1, 200):
        length = 7 * i
        print("length:", length)
        a.append(length)
        base = make_4_cadence_free(length)
        print(sum(base))
        b.append(sum(base))
        red = reduce_non_partial_4_cadences(base)
        print(sum(red))
        c.append(sum(red))
        red = reduce_non_partial_4_cadences(red)
        print(sum(red))
        d.append(sum(red))
        red = reduce_non_partial_4_cadences(red)
        print(sum(red))
        e.append(sum(red))

    plt.scatter(a, b)
    plt.show()

    plt.scatter(a, c)
    plt.show()

    plt.scatter(a, d)
    plt.show()

    plt.scatter(a, e)
    plt.show()


if __name__ == '__main__':
    test()