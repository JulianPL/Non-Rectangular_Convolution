#!/usr/bin/python3

import unittest

import nrconv

import sympy


class TestPowersOfTwo(unittest.TestCase):
    def test_create_power_of_two_simple1(self):
        result = nrconv.primes.create_power_of_two(6)
        want = 8
        self.assertEqual(result, want)

    def test_create_power_of_two_simple2(self):
        result = nrconv.primes.create_power_of_two(231)
        want = 256
        self.assertEqual(result, want)

    def test_create_power_of_two_power(self):
        result = nrconv.primes.create_power_of_two(64)
        want = 64
        self.assertEqual(result, want)

    def test_create_power_of_two_plus(self):
        result = nrconv.primes.create_power_of_two(33)
        want = 64
        self.assertEqual(result, want)

    def test_create_power_of_two_minus(self):
        result = nrconv.primes.create_power_of_two(15)
        want = 16
        self.assertEqual(result, want)

    def test_create_power_of_two_one(self):
        result = nrconv.create_power_of_two(1)
        want = 1
        self.assertEqual(result, want)


class TestPrimeGeneration(unittest.TestCase):
    def test_create_mod_prime_primality1(self):
        prime = nrconv.create_mod_prime(128, 1)
        self.assertTrue(sympy.ntheory.primetest.isprime(prime))

    def test_create_mod_prime_primality2(self):
        prime = nrconv.create_mod_prime(256, 1, 3121)
        self.assertTrue(sympy.ntheory.primetest.isprime(prime))

    def test_create_mod_prime_size(self):
        prime = nrconv.create_mod_prime(256, 1, 3121)
        self.assertGreaterEqual(prime, 3121)

    def test_create_mod_prime_modulo(self):
        mod = nrconv.create_mod_prime(256, 1, 3121) % 256
        remainder = 1
        self.assertEqual(mod, remainder)


class TestPrimesForNTT(unittest.TestCase):
    def test_create_ntt_prime_primality1(self):
        list1 = [1, 2, 3, 4, 5, 6, 7, 8]
        list2 = [8, 7, 6, 5, 4, 3, 2, 1]
        prime = nrconv.create_ntt_prime(list1, list2)
        self.assertTrue(sympy.ntheory.primetest.isprime(prime))

    def test_create_ntt_prime_primality2(self):
        list1 = [14, 23, 63, 41, 12, 42, 75, 32, 21]
        list2 = [14, 23, 63, 41, 12, 42, 75, 32, 21]
        prime = nrconv.create_ntt_prime(list1, list2)
        self.assertTrue(sympy.ntheory.primetest.isprime(prime))

    def test_create_ntt_prime_size(self):
        list1 = [14, 23, 63, 41, 12, 42, 75, 32, 21]
        list2 = [14, 23, 63, 41, 12, 42, 75, 32, 21]
        prime = nrconv.create_ntt_prime(list1, list2)
        self.assertGreaterEqual(prime, 15098 * 2 + 1)

    def test_create_ntt_prime_size(self):
        list1 = [42, 42, 42, 42, 42, 42, 42, 42, 42]
        list2 = [42, 42, 42, 42, 42, 42, 42, 42, 42]
        prime = nrconv.create_ntt_prime(list1, list2)
        self.assertGreaterEqual(prime, 15876 * 2 + 1)

    def test_create_ntt_prime_modulo(self):
        list1 = [14, 23, 63, 41, 12, 42, 75, 32]
        list2 = [14, 23, 63, 41, 12, 42, 75, 32]
        prime = nrconv.create_ntt_prime(list1, list2)
        mod = prime % 16
        remainder = 1
        self.assertEqual(mod, remainder)


if __name__ == '__main__':
    unittest.main()
