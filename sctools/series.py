from functools import lru_cache

from mpmath import beta, factorial, gamma, gammaprod, hyp2f1, nprod, nsum

from .tools import partition


def phi_series_coeffs(n: int, tau1, tau2, tau3) -> list:

    tau4 = 2 - (tau1 + tau2 + tau3)

    def a_coeff(n):
        return (
            (1 / factorial(n))
            * (gamma(tau1) / gamma(tau1 - n))
            * (beta(tau2 + n, tau3) / beta(tau2, tau3))
            * hyp2f1(n + tau2, 1 + n - tau1, n + tau2 + tau3, -1)
        )

    def b_coeff(n):
        return (
            (1 / factorial(n))
            * (gamma(tau1) / gamma(tau1 - n))
            * hyp2f1(tau4, 1 + n - tau1, tau3 + tau4, -1)
        )

    @lru_cache(maxsize=100)
    def c_coeff(n):
        if n == 0:
            return 1 / b_coeff(0)
        return -c_coeff(0) * nsum(lambda k: b_coeff(k) * c_coeff(n - k), [1, n])

    def d_coeff(n):
        return nsum(lambda k: a_coeff(n - k) * c_coeff(k), [0, n])

    k = beta(tau2, tau3) / beta(tau3, tau4)

    return [k * d_coeff(i) for i in range(n)]


def phi_series(x, coeffs):

    return sum(c * (x - 1) ** k for k, c in enumerate(coeffs))


def lam(*args):

    a = sum((j * k for j, k in enumerate(args, 2)))
    b = 1 + sum((j * k for j, k in enumerate(args, 1)))

    c = [k + 1 for k in args]

    return gammaprod([a + 1], [b + 1, *c])


def psi_coeff(n: int, phi_coeffs: list):

    if n == 0:
        return 1

    if n == 1:
        return 1 / phi_coeffs[1]

    n = n - 1

    parts = partition(n, n)

    def a(j):
        return phi_coeffs[j + 1] / phi_coeffs[1]

    return sum(
        (
            (-1) ** sum(k_args)
            * lam(*k_args)
            * nprod(lambda j: a(int(j)) ** k_args[int(j) - 1], [1, n])
            for k_args in parts
        )
    ) / phi_coeffs[1] ** (n + 1)


def psi_series(x, coeffs):

    return sum(c * (x - coeffs[0]) ** k for k, c in enumerate(coeffs))
