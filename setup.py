import setuptools
import re
import os
import sys

setuptools.setup(
    name='LARRY-dataset',
    version='0.0.1',
    python_requires='>3.7.0',
    author='Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard',
    author_email='mvinyard@broadinstitute.org',
    url = 'https://github.com/pinellolab/LARRY-dataset',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    description='LARRY Dataset: lineage and RNA recovery',
    packages = setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        "anndata>=0.8",
        "pytorch-lightning>=1.7.5",
        "wget>=3.2",
    ],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT"
)
