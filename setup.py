import setuptools
import re
import os
import sys

setuptools.setup(
    name='LARRY-dataset',
    version='0.0.2rc1',
    python_requires='>3.9.0',
    author='Michael E. Vinyard',
    author_email='mvinyard.ai@gmail.com',
    url = 'https://github.com/scDiffEq/LARRY-dataset',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    description='LARRY Dataset: lineage and RNA recovery',
    packages = setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        "anndata>=0.9.1",
  #      "pytorch-lightning>=2.0.2",
   #     "wget>=3.2",
    ],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT"
)
