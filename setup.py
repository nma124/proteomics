#!/usr/bin/env python3
"""
Setup script for Proteomics PRM Data Processing
"""

from setuptools import setup, find_packages
import pathlib

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="proteomics-prm-processing",
    version="1.0.0",
    description="Tools for processing Skyline PRM export data with heavy peptide dilution analysis",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Neev",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21.0", 
        "scikit-learn>=1.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "jupyter>=1.0.0",
        ],
        "docs": [
            "sphinx>=4.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "prm-process=main:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)