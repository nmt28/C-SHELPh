#!/usr/bin/env python
"""
Setup script for C-SHELPh.

Run using the following command:

$ pip install .

"""

from setuptools import setup
import glob

import cshelph

setup(
    name="cshelph",
    version=cshelph.CSHELPH_VERSION,
    description="C-SHELPh is the Classification of Sub-aquatic Height Extracted Photons. It is designed to isolate bathymetric photons in ICESat2 ATL03 files.",
    author="Nathan Thomas",
    author_email="Nathan.Thomas@edgehill.ac.uk",
    scripts=glob.glob("bin/*.py"),
    packages=["cshelph"],
    license="LICENSE.txt",
    url="https://github.com/nmt28/C-SHELPh",
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Remote Sensing Scientists",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
)
