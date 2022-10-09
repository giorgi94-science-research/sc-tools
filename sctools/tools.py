from collections import deque

from mpmath import beta, hyp2f1, sign, mpf


def phi(x, tau1, tau2, tau3):

    tau4 = 2 - (tau1 + tau2 + tau3)

    a = beta(tau2, tau3) * hyp2f1(tau3, 1 - tau1, tau2 + tau3, x / (1 + x))
    b = beta(tau3, tau4) * hyp2f1(tau3, 1 - tau1, tau3 + tau4, 1 / (1 + x))

    return a / b


def partition(n: int, k: int) -> list:

    if k == 1:
        return [(n,)]

    ans = []

    for j in range(n // k + 1):
        ans.extend([(*v, j) for v in partition(n - k * j, k - 1)])

    return ans


def find_root(f, err=1e-15, upper_bound=1000, digits=30):
    fl = 0
    fr = []

    def tonum(a, *b):
        d = "".join([str(i) for i in b])
        return mpf(str(a) + "." + d)

    for i in range(upper_bound + 1):
        if sign(f(fl)) != sign(f(fl + 1)):
            break
        else:
            fl += 1

    if i == upper_bound:
        return None

    if err is not None and abs(f(fl)) < err:
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


def calc_r_invariant(w1, w2, w3, w4, tau1, tau2, tau3, tau4):
    def __inner(args):
        (v1, t1), (v2, t2), (v3, t3), (v4, t4) = args

        phi0 = abs(v2 - v3) / abs(v3 - v4)

        r = find_root(lambda x: phi(x, t1, t2, t3) - phi0, upper_bound=2)

        return r

    data = deque([(w1, tau1), (w2, tau2), (w3, tau3), (w4, tau4)])

    for i in range(4):
        if data[0][1] > 1:
            data.rotate(-1)
            continue

        r = __inner(data)

        if r is not None:
            return r if i % 2 == 0 else 1 / r
        data.rotate(-1)

    return None
