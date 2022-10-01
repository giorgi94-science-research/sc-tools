from matplotlib import pyplot as plt
from sctools import Quadrilateral


def main():

    W = [-1 + 1j, 2 + 7j, 5 + 9j, 7 + 3j]

    vertices = [(w.real, w.imag) for w in W]

    X = [x for x, _ in vertices]
    Y = [y for _, y in vertices]
    X.append(X[0])
    Y.append(Y[0])

    quad = Quadrilateral(*vertices)

    r, A, C = quad.sc.calc()

    m = quad.sc.conformal_modulus(r)

    points = [quad.sc.quad(i - 5j, A, C, r) for i in range(-5, 5 + 1)]

    plt.scatter([p.real for p in points], [p.imag for p in points])

    points = [quad.sc.quad(i - 2j, A, C, r) for i in range(-5, 5 + 1)]

    plt.scatter([p.real for p in points], [p.imag for p in points])

    plt.plot(X, Y)
    plt.show()


main()
