"""
VPSC8 - Visco-Plastic Self-Consistent Crystal Plasticity Code

Python implementation of the VPSC8 polycrystal plasticity simulation code.
Supports multiple deformation mechanisms, hardening laws, and multi-phase materials.

Version: 8.0.0
Authors: Carlos N. Tomé, Ricardo A. Lebensohn, et al.
Python Port: 2024
License: See LICENSE file
"""

__version__ = "8.0.0"
__author__ = "Carlos N. Tomé, Ricardo A. Lebensohn, et al."

from .main import run_vpsc8

__all__ = ["run_vpsc8"]
