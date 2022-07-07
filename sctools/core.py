from mpmath import pi, log, hyp2f1, mpc, exp, beta, quad

from .tools import calc_r_invariant, phi


class SC:
    def __init__(
        self,
        w1: complex,
        w2: complex,
        w3: complex,
        w4: complex,
        tau1: float,
        tau2: float,
        tau3: float,
        tau4: float,
    ):
        self.w1, self.w2, self.w3, self.w4 = w1, w2, w3, w4
        self.tau1, self.tau2, self.tau3, self.tau4 = tau1, tau2, tau3, tau4

    @staticmethod
    def __hyp(a: float, b: float, c: float, z: float) -> float:

        return beta(a, b) * hyp2f1(a, 1 - c, a + b, z)

    def len_scale(self, r: float, theta: float):

        return (
            theta ** (-self.tau4)
            * (1 + theta) ** (1 - self.tau2)
            * (1 + theta + r * theta) ** (1 - self.tau3)
        )

    def len_1(self, r: float, theta: float = 1.0):
        c = self.len_scale(r, theta)

        return c * self.__hyp(self.tau4, self.tau1, self.tau3, -r)

    def len_2(self, r: float, theta: float = 1.0):
        c = self.len_scale(r, theta)
        return (
            r ** (self.tau3 - 1)
            * c
            * self.__hyp(self.tau2, self.tau1, self.tau3, -1 / r)
        )

    def len_3(self, r: float, theta: float = 1.0):
        c = self.len_scale(r, theta)
        return (
            r ** (self.tau2 + self.tau3 - 1)
            * c
            * self.__hyp(self.tau2, self.tau3, self.tau1, -r)
        )

    def len_4(self, r: float, theta: float = 1.0):
        c = self.len_scale(r, theta)
        return (
            r ** (-self.tau4) * c * self.__hyp(self.tau4, self.tau3, self.tau1, -1 / r)
        )

    def generalized_modulus(self, r: float):
        return phi(r, self.tau1, self.tau2, self.tau3)

    def conformal_modulus(self, r: float):
        return phi(r, 0.5, 0.5, 0.5)

    def calc(self):

        w1, w2, w3, w4 = self.w1, self.w2, self.w3, self.w4

        Cp = abs(w1 - w4) / (w1 - w4)
        Ap = -Cp * w4

        tau1, tau2, tau3, tau4 = self.tau1, self.tau2, self.tau3, self.tau4

        r = calc_r_invariant(w1, w2, w3, w4, tau1, tau2, tau3, tau4)

        if r is None:
            return None

        Cpp = (
            (2 ** (1 - tau2) / abs(w1 - w2))
            * (1 + 1 / (1 + r)) ** (1 - tau3)
            * beta(tau1, tau2)
            * hyp2f1(tau1, 1 - tau3, tau1 + tau2, 1 / (1 + r))
        )
        A = -Ap / Cp
        C = 1 / (Cp * Cpp)

        return r, A, C

    def get_vertices(self, A: complex, C: complex, r: float, theta: float = 1.0):
        def rotate(x):
            return exp(pi * 1j * x)

        tau1, tau2, tau3 = self.tau1, self.tau2, self.tau3

        s1 = 0
        s2 = rotate(tau1 - 1) * self.len_2(r, theta)
        s3 = rotate(tau1 + tau2 - 2) * self.len_3(r, theta)
        s4 = rotate(tau1 + tau2 + tau3 - 3) * self.len_4(r, theta)

        w1 = s1
        w2 = w1 + s2
        w3 = w2 + s3
        w4 = w3 + s4

        w1 -= w4
        w2 -= w4
        w3 -= w4
        w4 -= w4

        w1 = A + C * w1
        w2 = A + C * w2
        w3 = A + C * w3
        w4 = A + C * w4

        return [w1, w2, w3, w4]

    def error_calc(self, A: complex, C: complex, r: float, theta: float = 1.0):

        w1, w2, w3, w4 = self.w1, self.w2, self.w3, self.w4
        v1, v2, v3, v4 = self.get_vertices(A, C, r, theta)

        return abs(v1 - w1) + abs(v2 - w2) + abs(v3 - w3) + abs(v4 - w4)

    def quad(self, z: complex, A: complex, C: complex, r: float, theta: float = 1.0):

        tau1, tau2, tau3 = self.tau1, self.tau2, self.tau3

        return self.w1 + C * quad(
            lambda u: (1 - u) ** (tau1 - 1)
            * (1 - u / 2) ** (tau2 - 1)
            * (1 - u / (2 + r)) ** (tau3 - 1),
            [1, z],
        )


class Quadrilateral:
    def __init__(self, a, b, c, d):

        self.A, self.B, self.C, self.D = [mpc(a, b) for a, b in [a, b, c, d]]

        self.__set_tau()
        self.__set_lengths()

        self.sc = SC(
            self.A, self.B, self.C, self.D, self.tau1, self.tau2, self.tau3, self.tau4
        )

    @staticmethod
    def make_clockwise(vertices: list):
        # need to modify

        W = [mpc(a, b) for a, b in vertices]
        w1, w2, w3, w4 = W

        Cp = abs(w1 - w4) / (w1 - w4)
        Ap = -Cp * w4

        if (Ap + Cp * w2).imag > 0:
            vertices = vertices[::-1]

        return vertices

    @staticmethod
    def find_length(a: complex, b: complex):

        return abs(a - b)

    @staticmethod
    def find_angle(a: complex, b: complex, c: complex):

        v1 = a - b
        v2 = c - b

        _norm = abs(v1) * abs(v2)

        _cos = (v1.real * v2.real + v1.imag * v2.imag) / _norm
        _sin = (v1.real * v2.imag - v1.imag * v2.real) / _norm

        e = complex(_cos, _sin)

        angle = log(e).imag / pi

        return angle if angle > 0 else 2 + angle

    def __set_tau(self):
        w1, w2, w3, w4 = self.A, self.B, self.C, self.D

        angles = [
            self.find_angle(w4, w1, w2),
            self.find_angle(w1, w2, w3),
            self.find_angle(w2, w3, w4),
            self.find_angle(w3, w4, w1),
        ]

        assert abs(sum(angles) - 2) < 1e-8, "failed to find angles"

        self.tau1, self.tau2, self.tau3, self.tau4 = angles

    def __set_lengths(self):
        w1, w2, w3, w4 = self.A, self.B, self.C, self.D

        self.len1, self.len2, self.len3, self.len4 = (
            self.find_length(w4, w1),
            self.find_length(w1, w2),
            self.find_length(w2, w3),
            self.find_length(w3, w4),
        )

    def get_angles_degree(self):
        tau = [self.tau1, self.tau2, self.tau3, self.tau4]

        return [round(180 * a, 2) for a in tau]
