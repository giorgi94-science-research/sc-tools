# SC Tools

This package is based on our research to calculated Schwarz-Christoffel (SC) parameters and conformal modulus of a quadrilateral.

## Installation

Package can be installed using pip

```
pip install --user git+https://github.com/giorgi94/sc-tools.git
```

## Usage

To start using this package first you install it using `pip` or directly download it from github. To start working you can copy and paste following code

```python
from sctools import Quadrilateral, phi


"""
    After initializing Quadrilateral
    with clockwise oriented vertices
    it will automatically calculate
    its inner angles and side lengths.
"""
quad = Quadrilateral((1, 0), (-10, -15), (-8, 10), (15, 1))


# print its angles in degree
print(quad.get_angles_degree())

# print tau's, where pi*tau is inner angle
print(quad.tau1, quad.tau2, quad.tau3, quad.tau4)

# print side lengths
print(quad.len1, quad.len2, quad.len3, quad.len4)

"""
    find r, A, C parameters for the quadrilateral,
    which come from following mapping

    f(x) = A + C int_{-infty}^{x}
        (1-u)^(tau1-1) (1-u/2)^(tau2-1) (1-u/(2+r))^(tau3-1) du
"""
r, A, C = quad.sc.calc()

# print vertices calculated by SC mapping
print(quad.sc.get_vertices(A, C, r))

# print error between real and calculated vertices
print(quad.sc.error_calc(A, C, r))

# print conformal modulus of the quadrilateral
print(quad.sc.conformal_modulus(r))

```

### Generalized modulus

```python
from sctools import phi
import sctools.series as sc_series

# set some values for tau
tau1, tau2, tau3 = 1.27589191934767, 0.175999535986355, 0.406684994391987

# caclulate phi power series coefficients
phi_coeffs = sc_series.phi_series_coeffs(40, tau1, tau2, tau3)


"""
    following function prints real
    and approximate values for generalized modulus.
    Approximation is good when x is close to 1.
"""
def print_real_and_approx(x):

    print(phi(x, tau1, tau2, tau3))
    print(sc_series.phi_series(x, phi_coeffs))
    print()

print_real_and_approx(0.534)
print_real_and_approx(1.534)

print_real_and_approx(0.0534)
print_real_and_approx(3.534)


```

### Modular equation

A modular equation reduces to the generalized modular equation with signature $1/\tau$ and degree $p$ which can be rewritten in the form

$$
    \frac{{}_{2}F_{1}(1-\tau,\tau;1-\tau;1-s^2)}{{}_{2}F_{1}(1-\tau,\tau;1-\tau;s^2)}=p
    \frac{{}_{2}F_{1}(1-\tau,\tau;1-\tau;1-l^2)}{{}_{2}F_{1}(1-\tau,\tau;1-\tau;l^2)}
$$

which can be written using generalized modulus and we can solve it for $s$

```python
from mpmath import hyp2f1
import sctools.series as sc_series
from sctools import phi

tau = 0.43

phi_coeffs = sc_series.phi_series_coeffs(25, 1 - tau, tau, 1 - tau)
psi_coeffs = [sc_series.psi_coeff(i, phi_coeffs) for i in range(15)]

phi_tau = lambda x: phi(x, 1 - tau, tau, 1 - tau)
psi_tau = lambda y: sc_series.psi_series(y, psi_coeffs)

F = lambda x: hyp2f1(tau, 1 - tau, 1, 1 - x**2) / hyp2f1(tau, 1 - tau, 1, x**2)

l = 0.4
p = 0.5

s = (1 + psi_tau(p * phi_tau(1 / l**2 - 1))) ** (-0.5)

print(F(s))
print(p * F(l))

```