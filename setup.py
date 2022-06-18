from setuptools import setup
import re
import os
import sys

setup(
    name='LARRY-data',
    python_requires='>3.8.0',
    author='Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard',
    author_email='mvinyard@broadinstitute.org',
    url = 'https://github.com/pinellolab/LARRY-data',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    description='LARRY Dataset: lineage and RNA recovery',
    packages = [
        'larry',
    ],
    install_requires=[
        "anndata>=0.7.1",
        "scanpy>=1.4.3",
        "torch>=1.1.0",
        "numpy>=1.19.2",
        "pandas>=1.1.2",
    ],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT"
)
