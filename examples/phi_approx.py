import sctools.series as sc_series
from sctools import phi


def phi_rect_approx(x):
    return (
        1
        + 0.2284732905 * (x - 1)
        - 0.08813662302 * (x - 1) ** 2
        + 0.04827946003 * (x - 1) ** 3
        - 0.03127277573 * (x - 1) ** 4
        + 0.02241699206 * (x - 1) ** 5
        - 0.01718041161 * (x - 1) ** 6
        + 0.01379141315 * (x - 1) ** 7
        - 0.01144586284 * (x - 1) ** 8
        + 0.009737374309 * (x - 1) ** 9
        - 0.008442674960 * (x - 1) ** 10
        + 0.007430490326 * (x - 1) ** 11
        - 0.006619174911 * (x - 1) ** 12
        + 0.005955504769 * (x - 1) ** 13
        - 0.005403384294 * (x - 1) ** 14
        + 0.004937503063 * (x - 1) ** 15
    )


def inv_phi_rect_approx(x):

    return (
        1.0
        + 4.376879230452952 * (x - 1)
        + 7.390096283758735 * (x - 1) ** 2
        + 7.237219471534163 * (x - 1) ** 3
        + 5.991576168110211 * (x - 1) ** 4
        + 3.728577509035313 * (x - 1) ** 5
        + 1.7396270463685406 * (x - 1) ** 6
        + 1.1543413391162676 * (x - 1) ** 7
        + 0.04250082824813188 * (x - 1) ** 8
        + 0.3292792774810219 * (x - 1) ** 9
        - 0.002590010314452334 * (x - 1) ** 10
        - 0.16126267764293367 * (x - 1) ** 11
        + 0.36744094492367263 * (x - 1) ** 12
        - 0.4756371436446193 * (x - 1) ** 13
        + 0.4626527213655051 * (x - 1) ** 14
        - 0.3114696848504541 * (x - 1) ** 15
    )


def series_test():

    tau = 0.75

    tau1, tau2, tau3 = 1 - tau, tau, 1 - tau

    tau4 = 2 - (tau1 + tau2 + tau3)

    assert tau4 > 0, "angles are incorrect"

    phi_n_list = sc_series.phi_series_coeffs(10, tau1, tau2, tau3)

    x = 1.34

    print("approx\t", sc_series.phi_series(x, phi_n_list))
    print("exact\t", phi(x, tau1, tau2, tau3))

    y = phi(x, tau1, tau2, tau3)

    inv_phi_n_list = [sc_series.psi_coeff(i, phi_n_list) for i in range(5)]

    print("inv\t", sc_series.phi_series(y, inv_phi_n_list))
