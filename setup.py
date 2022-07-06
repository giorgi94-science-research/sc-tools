from pathlib import Path

from setuptools import setup

here = Path(__file__).resolve().parent


README = (here / "README.md").read_text(encoding="utf-8")

requirements = [
    line.strip()
    for line in (here / "requirements.txt").read_text(encoding="utf-8").split("\n")
    if line.strip()
]


setup(
    name="sc-tools",
    version="2.0",
    packages=["sctools"],
    description="Tools to work with Schwarz-Christoffel mapping",
    include_package_data=True,
    long_description=README,
    long_description_content_type="text/markdown",
    author="Giorgi Kakulashvili",
    url="https://github.com/giorgi94/sc-tools",
    keywords=["sc", "Schwarz-Christoffel", "conformal modulus", "quadrilateral"],
    platforms=["OS Independent"],
    install_requires=requirements,
)
