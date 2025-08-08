"""
tensor.py
---------
Tensor algebra utilities for VPSC Python port.
Implements Voigt notation conversions, basis changes, and matrix inversion.
"""

import numpy as np


def voigt(T1=None, T2=None, C2=None, C4=None, opt=1):
    """
    Implements the VOIGT subroutine from VPSC FORTRAN.
    Converts between Voigt and tensor notation for 2nd and 4th order tensors.
    Args:
        T1: 6-vector (Voigt)
        T2: 3x3 tensor
        C2: 6x6 matrix (Voigt)
        C4: 3x3x3x3 tensor
        opt: operation mode (1,2,3,4)
    Returns:
        Tuple of converted tensors as needed.
    """
    ijv = np.array([[0, 0], [1, 1], [2, 2], [1, 2], [0, 2], [0, 1]])
    if opt == 1:
        # T1 (6,) -> T2 (3,3)
        T2 = np.zeros((3, 3))
        for i in range(6):
            I1, I2 = ijv[i]
            T2[I1, I2] = T1[i]
            T2[I2, I1] = T1[i]
        return T2
    elif opt == 2:
        # T2 (3,3) -> T1 (6,)
        T1 = np.zeros(6)
        for i in range(6):
            I1, I2 = ijv[i]
            T1[i] = T2[I1, I2]
        return T1
    elif opt == 3:
        # C2 (6,6) -> C4 (3,3,3,3)
        C4 = np.zeros((3, 3, 3, 3))
        for i in range(6):
            I1, I2 = ijv[i]
            for j in range(6):
                J1, J2 = ijv[j]
                C4[I1, I2, J1, J2] = C2[i, j]
                C4[I2, I1, J1, J2] = C2[i, j]
                C4[I1, I2, J2, J1] = C2[i, j]
                C4[I2, I1, J2, J1] = C2[i, j]
        return C4
    elif opt == 4:
        # C4 (3,3,3,3) -> C2 (6,6)
        C2 = np.zeros((6, 6))
        for i in range(6):
            I1, I2 = ijv[i]
            for j in range(6):
                J1, J2 = ijv[j]
                C2[i, j] = C4[I1, I2, J1, J2]
        return C2
    else:
        raise ValueError("Invalid opt for voigt")


def lu_inverse(A):
    """
    Inverts a matrix using LU decomposition (NumPy wrapper).
    Normalizes for numerical stability as in FORTRAN code.
    Args:
        A: (n,n) array
    Returns:
        Ainv: (n,n) array, inverse of A
    Raises:
        np.linalg.LinAlgError if singular.
    """
    A = np.array(A, dtype=float)
    amax = np.max(np.abs(A))
    if amax == 0:
        raise np.linalg.LinAlgError("Singular matrix in lu_inverse")
    A_norm = A / amax
    try:
        Ainv = np.linalg.inv(A_norm)
    except np.linalg.LinAlgError:
        raise np.linalg.LinAlgError("Singular matrix in lu_inverse")
    return Ainv / amax


def chg_basis(CE2=None, C2=None, CE4=None, C4=None, opt=None, kdim=6):
    """
    Implements the CHG_BASIS subroutine from VPSC FORTRAN.
    Converts between b-basis (Voigt-like) and full tensor notation for 2nd and 4th order tensors.
    Args:
        CE2: (kdim,) vector in b-basis
        C2: (3,3) tensor
        CE4: (kdim,kdim) matrix in b-basis
        C4: (3,3,3,3) tensor
        opt: operation mode (0-4)
        kdim: basis dimension (default 6)
    Returns:
        Tuple of converted tensors as needed.
    """
    import math

    B = np.zeros((3, 3, kdim))
    sqrt2 = math.sqrt(2.0)
    sqrt3 = math.sqrt(3.0)
    sqrt6 = math.sqrt(6.0)

    # Option 0: Initialize B tensor
    if opt == 0:
        B[:, :, 1] = 0
        B[:, :, 2] = 0
        B[:, :, 3] = 0
        B[:, :, 4] = 0
        B[:, :, 5] = 0
        B[:, :, 0] = 0
        B[0, 0, 1] = -1.0 / sqrt6
        B[1, 1, 1] = -1.0 / sqrt6
        B[2, 2, 1] = 2.0 / sqrt6
        B[0, 0, 0] = -1.0 / sqrt2
        B[1, 1, 0] = 1.0 / sqrt2
        B[1, 2, 2] = 1.0 / sqrt2
        B[2, 1, 2] = 1.0 / sqrt2
        B[0, 2, 3] = 1.0 / sqrt2
        B[2, 0, 3] = 1.0 / sqrt2
        B[0, 1, 4] = 1.0 / sqrt2
        B[1, 0, 4] = 1.0 / sqrt2
        B[0, 0, 5] = 1.0 / sqrt3
        B[1, 1, 5] = 1.0 / sqrt3
        B[2, 2, 5] = 1.0 / sqrt3
        return B

    # Always initialize B for other options
    B = np.zeros((3, 3, kdim))
    B[0, 0, 1] = -1.0 / sqrt6
    B[1, 1, 1] = -1.0 / sqrt6
    B[2, 2, 1] = 2.0 / sqrt6
    B[0, 0, 0] = -1.0 / sqrt2
    B[1, 1, 0] = 1.0 / sqrt2
    B[1, 2, 2] = 1.0 / sqrt2
    B[2, 1, 2] = 1.0 / sqrt2
    B[0, 2, 3] = 1.0 / sqrt2
    B[2, 0, 3] = 1.0 / sqrt2
    B[0, 1, 4] = 1.0 / sqrt2
    B[1, 0, 4] = 1.0 / sqrt2
    B[0, 0, 5] = 1.0 / sqrt3
    B[1, 1, 5] = 1.0 / sqrt3
    B[2, 2, 5] = 1.0 / sqrt3

    # Option 1: CE2 (kdim,) -> C2 (3,3)
    if opt == 1:
        if CE2 is None:
            raise ValueError("CE2 must be provided for opt=1")
        C2 = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                for n in range(kdim):
                    C2[i, j] += CE2[n] * B[i, j, n]
        return C2

    # Option 2: C2 (3,3) -> CE2 (kdim,)
    if opt == 2:
        if C2 is None:
            raise ValueError("C2 must be provided for opt=2")
        CE2 = np.zeros(kdim)
        for n in range(kdim):
            for i in range(3):
                for j in range(3):
                    CE2[n] += C2[i, j] * B[i, j, n]
        return CE2

    # Option 3: CE4 (kdim,kdim) -> C4 (3,3,3,3)
    if opt == 3:
        if CE4 is None:
            raise ValueError("CE4 must be provided for opt=3")
        C4 = np.zeros((3, 3, 3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for n in range(kdim):
                            for m in range(kdim):
                                C4[i, j, k, l] += CE4[n, m] * B[i, j, n] * B[k, l, m]
        return C4

    # Option 4: C4 (3,3,3,3) -> CE4 (kdim,kdim)
    if opt == 4:
        if C4 is None:
            raise ValueError("C4 must be provided for opt=4")
        CE4 = np.zeros((kdim, kdim))
        for n in range(kdim):
            for m in range(kdim):
                for i in range(3):
                    for j in range(3):
                        for k in range(3):
                            for l in range(3):
                                CE4[n, m] += C4[i, j, k, l] * B[i, j, n] * B[k, l, m]
        return CE4

    raise ValueError("Invalid opt for chg_basis")
