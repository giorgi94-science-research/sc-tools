import operator
from functools import lru_cache, reduce
from math import ceil

from mpmath import beta, gamma, hyper, mp, quad, pi, log, acos, sign

mp.dps = 20


def prod(a, start=1):
    return reduce(operator.mul, a, start)


def factorial(n):
    return gamma(n + 1)


def integrate(f, a, b):
    return quad(f, [a, b])


def make_clockwise(vertices: list):

    W = [a + 1j * b for a, b in vertices]
    w1, w2, w3, w4 = W

    Cp = abs(w1 - w4) / (w1 - w4)
    Ap = -Cp * w4

    if (Ap + Cp * w2).imag > 0:
        vertices = vertices[::-1]

    return vertices


def find_angle(a: complex, b: complex, c: complex):
    v1 = a - b
    v2 = c - b

    _norm = abs(v1) * abs(v2)

    _cos = (v1.real * v2.real + v1.imag * v2.imag) / _norm
    _sin = (v1.real * v2.imag - v1.imag * v2.real) / _norm

    alpha = acos(abs(_cos)) / pi

    alpha = (log(_cos + 1j * _sin) / (pi * 1j)).real

    if alpha < 0:
        alpha = 2 + alpha

    return alpha


def find_quad_angles(vertices: list):

    W = [a + 1j * b for a, b in vertices]
    w1, w2, w3, w4 = W

    angles = [
        find_angle(w4, w1, w2),
        find_angle(w1, w2, w3),
        find_angle(w2, w3, w4),
        find_angle(w3, w4, w1),
    ]

    assert abs(sum(angles) - 2) < 1e-8, "failed to find angles"

    return angles


def find_root(f, err=1e-15, upper_bound=10, digits=30):
    fl = 0
    fr = []

    def tonum(a, *b):
        d = "".join([str(i) for i in b])
        return eval(str(a) + "." + d)

    for _ in range(ceil(upper_bound)):
        if sign(f(fl)) != sign(f(fl + 1)):
            break
        else:
            fl += 1

    if abs(f(fl)) < err:
        return fl

    for _ in range(digits):
        d = 0
        for _ in range(9):
            a = tonum(fl, *fr, d)
            b = tonum(fl, *fr, d + 1)

            if sign(f(a)) != sign(f(b)):
                fr.append(d)
                break
            else:
                d += 1
        if d == 9:
            fr.append(d)

    return tonum(fl, *fr)


def F(a: float, b: float, c: float, z: float) -> float:
    return hyper([a, b], [c], z)


def G(a: float, b: float, c: float, z: float) -> float:

    return beta(a, b) * F(a, 1 - c, a + b, z)


def side_length_scale(r: float, theta: float, tau1: float, tau2: float, tau3: float):
    tau4 = 2 - (tau1 + tau2 + tau3)

    return (
        theta ** (-tau4)
        * (1 + theta) ** (1 - tau2)
        * (1 + theta + r * theta) ** (1 - tau3)
    )


def side_length_1(r: float, theta: float, tau1: float, tau2: float, tau3: float):
    tau4 = 2 - (tau1 + tau2 + tau3)
    c = side_length_scale(r, theta, tau1, tau2, tau3)

    return c * G(tau4, tau1, tau3, -r)


def side_length_2(r: float, theta: float, tau1: float, tau2: float, tau3: float):
    c = side_length_scale(r, theta, tau1, tau2, tau3)
    return r ** (tau3 - 1) * c * G(tau2, tau1, tau3, -1 / r)


def side_length_3(r: float, theta: float, tau1: float, tau2: float, tau3: float):
    c = side_length_scale(r, theta, tau1, tau2, tau3)
    return r ** (tau2 + tau3 - 1) * c * G(tau2, tau3, tau1, -r)


def side_length_4(r: float, theta: float, tau1: float, tau2: float, tau3: float):
    tau4 = 2 - (tau1 + tau2 + tau3)
    c = side_length_scale(r, theta, tau1, tau2, tau3)
    return r ** (-tau4) * c * G(tau4, tau3, tau1, -1 / r)


def phi(r, tau1, tau2, tau3):
    if r == 0:
        return 0

    tau4 = 2 - (tau1 + tau2 + tau3)

    return (
        r ** (1 - tau1)
        * (beta(tau2, tau3) / beta(tau3, tau4))
        * F(tau2, 1 - tau1, tau2 + tau3, -r)
        / F(tau4, 1 - tau1, tau3 + tau4, -1 / r)
    )


def K(r):
    return (pi / 2) * F(0.5, 0.5, 1, r ** 2)


def conformal_moduli(r):
    k = (1 / (1 + r)) ** 0.5
    return K((1 - k ** 2) ** 0.5) / K(k)


def calc_sc_params(vertices, tau, upper_bound=50):
    W = [a + 1j * b for a, b in vertices]

    w1, w2, w3, w4 = W

    Cp = abs(w1 - w4) / (w1 - w4)
    Ap = -Cp * w4

    tau1, tau2, tau3, _ = tau

    phi0 = abs(w2 - w3) / abs(w3 - w4)

    def _f(x):
        return phi(x, tau1, tau2, tau3) - phi0

    r = find_root(lambda x: _f(x), upper_bound=upper_bound)

    Cpp = (
        (2 ** (1 - tau2) / abs(w1 - w2))
        * (1 + 1 / (1 + r)) ** (1 - tau3)
        * beta(tau1, tau2)
        * F(tau1, 1 - tau3, tau1 + tau2, 1 / (1 + r))
    )
    A = -Ap / Cp
    C = 1 / (Cp * Cpp)

    return r, A, C


def phi_series_coeffs(n: int, tau1, tau2, tau3) -> list:

    tau4 = 2 - (tau1 + tau2 + tau3)

    def a_coeff(n):
        return (
            (1 / factorial(n))
            * (gamma(tau1) / gamma(tau1 - n))
            * (beta(tau2 + n, tau3) / beta(tau2, tau3))
            * F(n + tau2, 1 + n - tau1, n + tau2 + tau3, -1)
        )

    def b_coeff(n):
        return (
            (1 / factorial(n))
            * (gamma(tau1) / gamma(tau1 - n))
            * F(tau4, 1 + n - tau1, tau3 + tau4, -1)
        )

    @lru_cache(maxsize=100)
    def c_coeff(n):
        if n == 0:
            return 1 / b_coeff(0)
        return -c_coeff(0) * sum(b_coeff(k) * c_coeff(n - k) for k in range(1, n + 1))

    def d_coeff(n):
        return sum(a_coeff(n - k) * c_coeff(k) for k in range(0, n + 1))

    k = beta(tau2, tau3) / beta(tau3, tau4)

    return [k * d_coeff(i) for i in range(n)]


def phi_series(x, coeffs):

    return sum(c * (x - 1) ** k for k, c in enumerate(coeffs))


def partition(n: int, k: int) -> list:
    result = []

    def step(*args):
        j = len(args)

        s = sum((k - i) * a for i, a in enumerate(reversed(args)))

        if j == k - 1:
            return result.append((n - s, *args))

        m = (n - s) // (k - j)

        for a in range(m + 1):
            step(a, *args)

    k > 0 and step()

    return result


def psi_series_coeffs(phi_coeffs: list):
    def lam(*k):
        _a = factorial(sum((i + 2) * ki for i, ki in enumerate(k))) / factorial(
            1 + sum((i + 1) * ki for i, ki in enumerate(k))
        )
        for i in k:
            _a /= factorial(i)
        return _a

    def c(n):
        if n == 0:
            return phi_coeffs[0]
        if n == 1:
            return 1 / phi_coeffs[1]

        return (
            sum(
                (-1) ** sum(k)
                * lam(*k)
                * prod(
                    (phi_coeffs[j + 2] / phi_coeffs[1]) ** kj for j, kj in enumerate(k)
                )
                for k in partition(n - 1, n - 1)
            )
            / phi_coeffs[1] ** n
        )

    return [c(n) for n in range(len(phi_coeffs) - 1)]


def psi_series(x, coeffs):

    return sum(c * (x - coeffs[0]) ** k for k, c in enumerate(coeffs))
