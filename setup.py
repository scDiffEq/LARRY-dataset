# -- import packages: ---------------------------------------------------------
import setuptools
import re
import os
import sys

# -- fetch requirements packages: ---------------------------------------------
with open("requirements.txt") as f:
    requirements = f.read().splitlines()

with open(f"larry/__version__.py") as v:
    exec(v.read())

# -- run setup: ---------------------------------------------------------------
setuptools.setup(
    name='LARRY-dataset',
    version='0.0.2rc1',
    python_requires='>3.10',
    author='Michael E. Vinyard',
    author_email='mvinyard.ai@gmail.com',
    url = 'https://github.com/scDiffEq/LARRY-dataset',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    description='LARRY Dataset: lineage and RNA recovery',
    packages = setuptools.find_packages(),
    include_package_data=True,
    install_requires=requirements,
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.10",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT"
)
