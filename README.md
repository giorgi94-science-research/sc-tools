# SC Tools

This package is based on my research to calculated Schwarz-Christoffel (SC) parameters and conformal modulus of a quadrilateral.

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