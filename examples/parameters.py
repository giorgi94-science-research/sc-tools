from mpmath import hyp2f1
from sctools import Quadrilateral, phi

import sctools.series as sc_series


def compare_conformal_moduli_1():
    """
    compare with https://arxiv.org/pdf/2103.10237.pdf
    page 16, table 1
    """

    compare = [
        [7 + 5j, -1 + 2j, 1.17336589158553],
        [8 + 3j, -1 + 1j, 0.71853428024898],
        [5 + 5j, -3 + 1j, 1.00171178298845],
        [7 + 4j, -3 + 3j, 1.17821610141750],
        [5 + 5j, -1 + 2j, 1.27382477147819],
        [7 + 5j, 0 + 1j, 0.92223220304256],
        [7 + 3j, 1 + 2j, 1.68574560877551],
        [4 + 5j, -2 + 1j, 1.02479880902234],
    ]

    for A, B, module in compare:
        vertices = [(B.real, B.imag), (A.real, A.imag), (1, 0), (0, 0)]

        quad = Quadrilateral(*vertices)

        r, A, C = quad.sc.calc()

        m = quad.sc.conformal_modulus(r)

        print("caclulated module:", "\t", m)
        print("compare module:", "\t", module)

        print()


def compare_conformal_moduli_2():
    """
    Heikkala, V., Vamanamurthy, M.K. & Vuorinen,
    M. Generalized Elliptic Integrals.
    Comput. Methods Funct. Theory 9, 75–109 (2009).
    page 84, table 1
    """

    def check(m, n):
        vertices = [
            1 + 0j,
            0j,
            1j,
            m + 1j * n,
        ]

        points = [(v.real, v.imag) for v in vertices]

        quad = Quadrilateral(*points)

        r, A, C = quad.sc.calc()

        return quad.sc.conformal_modulus(r)

    for m in range(1, 6):
        for n in range(1, 6):
            print(f"m={m}, n={n} :", check(m, n))
        print()


def modular_equation():
    a = 0.43

    phi_coeffs = sc_series.phi_series_coeffs(25, 1 - a, a, 1 - a)
    psi_coeffs = [sc_series.psi_coeff(i, phi_coeffs) for i in range(15)]

    phi_a = lambda x: phi(x, 1 - a, a, 1 - a)
    psi_a = lambda y: sc_series.psi_series(y, psi_coeffs)

    F = lambda x: hyp2f1(a, 1 - a, 1, 1 - x**2) / hyp2f1(a, 1 - a, 1, x**2)

    r = 0.4
    p = 0.5

    s = (1 + psi_a(p * phi_a(1 / r**2 - 1))) ** (-0.5)

    print(F(s))
    print(p * F(r))


if __name__ == "__main__":

    # compare_conformal_moduli_1()
    # compare_conformal_moduli_2()

    modular_equation()
