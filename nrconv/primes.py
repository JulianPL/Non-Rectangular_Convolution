#!/usr/bin/python3
"""A module to create suitable primes for the NTT.
"""

from typing import List

import sympy


def create_power_of_two(number: int) -> int:
    """Creates the smallest power of two which is at least number."""
    power = 1
    while power < number:
        power = power * 2
    return power


def create_mod_prime(base: int, mod: int, min_prime: int = 0) -> int:
    """Creates a prime >= min_prime with prime % base == mod

    The algorithm is deterministic.  In the worst case, at most
    ((max*log(max))^2)/len(list) numbers need to be tested.

    This algorithm might be improvable with an adapted sieve
    of Eratosthenes:  For details see the appendix of:
    https://drops.dagstuhl.de/opus/volltexte/2020/11891/pdf/LIPIcs-STACS-2020-30.pdf
    Also, random numbers might be helpful.  For details see section 2.3 of:
    https://epubs.siam.org/doi/pdf/10.1137/100811167
    """

    prime_candidate = base + mod
    if prime_candidate < min_prime:
        diff_number_bases = (min_prime - prime_candidate + base - 1) // base
        prime_candidate = prime_candidate + diff_number_bases * base
    while not sympy.ntheory.primetest.isprime(prime_candidate):
        prime_candidate = prime_candidate + base
    return prime_candidate


def create_ntt_prime(list1: List[int], list2: List[int]) -> int:
    """Creates a suitable prime for the NTT.

    The prime is of the form m*2^k + 1, where 2^k is the length of the NTT.
    Furthermore, the prime is guaranteed to be big enough in order to
    prevent overflows in the NTT.
    """

    ntt_length = 2 * create_power_of_two(max(len(list1), len(list2)))
    max_list1_abs = max(abs(max(list1)), abs(min(list1)))
    max_list2_abs = max(abs(max(list2)), abs(min(list2)))
    max_value = max_list1_abs * max_list2_abs * ntt_length + 1

    prime = create_mod_prime(ntt_length, 1, max_value)
    return prime
