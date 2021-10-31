from mpmath import hyper
from sctools.core import (
    calc_sc_params,
    conformal_moduli,
    find_quad_angles,
    make_clockwise,
    phi,
    phi_series_coeffs,
    psi_series,
    psi_series_coeffs,
    side_length_1,
    side_length_2,
    side_length_3,
    side_length_4,
)


def calc_conformal_moduli(vertices):

    vertices = make_clockwise(vertices)

    tau = find_quad_angles(vertices)

    r, A, C = calc_sc_params(vertices, tau, upper_bound=1000)

    return conformal_moduli(r)


def compare_conformal_moduli():
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
        vertices = [(0, 0), (1, 0), (A.real, A.imag), (B.real, B.imag)]

        m = calc_conformal_moduli(vertices)

        print("my module:", m)
        print("co module:", module)

        print()


def check_side_lengths(m=2, n=4):

    W = [m + 1j * n, 1j, 0j, 1 + 0j]

    vertices = [(z.real, z.imag) for z in W]
    vertices = make_clockwise(vertices)
    W = [a + 1j * b for a, b in vertices]

    tau1, tau2, tau3, tau4 = find_quad_angles(vertices)

    r, A, C = calc_sc_params(vertices, (tau1, tau2, tau3, tau4), upper_bound=50)

    s1 = abs(W[3] - W[0])
    s2 = abs(W[0] - W[1])
    s3 = abs(W[1] - W[2])
    s4 = abs(W[2] - W[3])

    scale = abs(C)

    print(s1, scale * side_length_1(r, 1, tau1, tau2, tau3))
    print(s2, scale * side_length_2(r, 1, tau1, tau2, tau3))
    print(s3, scale * side_length_3(r, 1, tau1, tau2, tau3))
    print(s4, scale * side_length_4(r, 1, tau1, tau2, tau3))


def modular_equation():
    a = 0.43

    phi_coeffs = phi_series_coeffs(25, 1 - a, a, 1 - a)
    psi_coeffs = psi_series_coeffs(phi_coeffs)

    phi_a = lambda x: phi(x, 1 - a, a, 1 - a)
    psi_a = lambda y: psi_series(y, psi_coeffs)

    F = lambda x: hyper([a, 1 - a], [1], 1 - x ** 2) / hyper([a, 1 - a], [1], x ** 2)

    r = 0.4
    p = 0.5

    s = (1 + psi_a(p * phi_a(1 / r ** 2 - 1))) ** (-0.5)

    print(F(s))
    print(p * F(r))
