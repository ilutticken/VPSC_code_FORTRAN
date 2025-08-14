"""
Setup script for VPSC8 - Visco-Plastic Self-Consistent Crystal Plasticity Code
"""

from setuptools import setup, find_packages
import os


# Read the README file
def read_readme():
    try:
        with open("README.md", "r", encoding="utf-8") as fh:
            return fh.read()
    except FileNotFoundError:
        return "VPSC8 - Visco-Plastic Self-Consistent Crystal Plasticity Code"


# Read requirements
def read_requirements():
    try:
        with open("requirements.txt", "r", encoding="utf-8") as fh:
            return [
                line.strip() for line in fh if line.strip() and not line.startswith("#")
            ]
    except FileNotFoundError:
        return ["numpy>=1.19.0", "scipy>=1.5.0"]


setup(
    name="vpsc8",
    version="8.0.0",
    author="Carlos N. Tomé, Ricardo A. Lebensohn, et al.",
    author_email="",
    description="Visco-Plastic Self-Consistent Crystal Plasticity Code",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/lanl/VPSC8",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Materials Science",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
        "plotting": [
            "matplotlib>=3.3.0",
            "plotly>=5.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "vpsc8=vpsc8.main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "vpsc8": ["*.txt", "*.dat"],
    },
    zip_safe=False,
)
