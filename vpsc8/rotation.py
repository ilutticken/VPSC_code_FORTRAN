"""
Texture rotation functionality for VPSC8.
"""

import numpy as np


def read_rotation_matrix(filename):
    """
    Read a 3x3 rotation matrix from a rotation file.

    Parameters:
    -----------
    filename : str
        Path to the rotation matrix file

    Returns:
    --------
    numpy.ndarray
        3x3 rotation matrix
    """
    import os

    # Handle relative paths by looking in parent directory if needed
    original_filename = filename
    if not os.path.exists(filename):
        # Try parent directory
        parent_path = os.path.join("..", filename)
        if os.path.exists(parent_path):
            filename = parent_path
        else:
            # Try extracting just the filename and looking in parent
            basename = os.path.basename(filename)
            test_paths = [
                os.path.join("..", basename),
                os.path.join("..", "..", basename),
                os.path.join("..", "..", "..", basename),
                os.path.join(
                    "..", "..", filename.replace("\\", os.sep).replace("/", os.sep)
                ),
                os.path.join(
                    "..",
                    "..",
                    "..",
                    filename.replace("\\", os.sep).replace("/", os.sep),
                ),
            ]
            for test_path in test_paths:
                if os.path.exists(test_path):
                    filename = test_path
                    break

    try:
        with open(filename, "r") as f:
            lines = f.readlines()

        # Find the rotation matrix (skip comment lines)
        matrix_lines = []
        for line in lines:
            line = line.strip()
            if not line.startswith("*") and not line.startswith("c") and line:
                tokens = line.split()
                if len(tokens) >= 3:
                    try:
                        row = [float(tokens[0]), float(tokens[1]), float(tokens[2])]
                        matrix_lines.append(row)
                    except ValueError:
                        continue

        if len(matrix_lines) != 3:
            raise ValueError(
                f"Expected 3x3 rotation matrix, got {len(matrix_lines)} rows"
            )

        rotation_matrix = np.array(matrix_lines)

        # Validate that it's approximately a rotation matrix
        # Check determinant is close to 1
        det = np.linalg.det(rotation_matrix)
        if abs(det - 1.0) > 1e-6:
            print(f"Warning: Rotation matrix determinant is {det}, expected 1.0")

        # Check orthogonality: R^T * R should be identity
        product = rotation_matrix.T @ rotation_matrix
        identity = np.eye(3)
        if not np.allclose(product, identity, atol=1e-6):
            print("Warning: Rotation matrix is not orthogonal")

        return rotation_matrix

    except Exception as e:
        raise RuntimeError(
            f"Error reading rotation matrix from {original_filename}: {e}"
        )


def apply_texture_rotation(config, rotation_matrix):
    """
    Apply a rotation matrix to the texture (grain orientations).

    Parameters:
    -----------
    config : dict
        VPSC configuration containing grain orientations
    rotation_matrix : numpy.ndarray
        3x3 rotation matrix to apply
    """
    try:
        # Apply rotation to each grain's orientation
        for i_phase in range(config.get("nph", 1)):
            phase_key = f"phase_{i_phase+1}"
            if phase_key in config:
                phase_data = config[phase_key]

                if "orientations" in phase_data:
                    orientations = phase_data["orientations"]
                    n_grains = len(orientations)

                    print(
                        f"  Applying rotation to {n_grains} grains in phase {i_phase+1}"
                    )

                    for i_grain in range(n_grains):
                        # Get current orientation matrix (grain to sample transformation)
                        euler_angles = orientations[i_grain]

                        # Convert Euler angles to rotation matrix
                        phi1, theta, phi2 = np.radians(euler_angles)

                        # Bunge convention: Z-X-Z rotation
                        cos_phi1, sin_phi1 = np.cos(phi1), np.sin(phi1)
                        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
                        cos_phi2, sin_phi2 = np.cos(phi2), np.sin(phi2)

                        # Original orientation matrix
                        g_matrix = np.array(
                            [
                                [
                                    cos_phi1 * cos_phi2
                                    - sin_phi1 * sin_phi2 * cos_theta,
                                    sin_phi1 * cos_phi2
                                    + cos_phi1 * sin_phi2 * cos_theta,
                                    sin_phi2 * sin_theta,
                                ],
                                [
                                    -cos_phi1 * sin_phi2
                                    - sin_phi1 * cos_phi2 * cos_theta,
                                    -sin_phi1 * sin_phi2
                                    + cos_phi1 * cos_phi2 * cos_theta,
                                    cos_phi2 * sin_theta,
                                ],
                                [
                                    sin_phi1 * sin_theta,
                                    -cos_phi1 * sin_theta,
                                    cos_theta,
                                ],
                            ]
                        )

                        # Apply rotation: g_new = R * g_old
                        g_rotated = rotation_matrix @ g_matrix

                        # Convert back to Euler angles
                        new_euler = rotation_matrix_to_euler(g_rotated)

                        # Update the orientation
                        orientations[i_grain] = np.degrees(new_euler)

    except Exception as e:
        print(f"Error applying texture rotation: {e}")
        raise


def rotation_matrix_to_euler(R):
    """
    Convert a 3x3 rotation matrix to Euler angles (Bunge convention: Z-X-Z).

    Parameters:
    -----------
    R : numpy.ndarray
        3x3 rotation matrix

    Returns:
    --------
    numpy.ndarray
        Euler angles [phi1, theta, phi2] in radians
    """
    # Extract Euler angles from rotation matrix (Bunge convention)
    theta = np.arccos(np.clip(R[2, 2], -1.0, 1.0))

    if abs(R[2, 2]) >= 0.9999:  # theta ≈ 0 or π
        phi2 = 0.0
        phi1 = np.arctan2(R[0, 1], R[0, 0])
    else:
        sin_theta = np.sin(theta)
        phi2 = np.arctan2(R[0, 2] / sin_theta, R[1, 2] / sin_theta)
        phi1 = np.arctan2(R[2, 0] / sin_theta, -R[2, 1] / sin_theta)

    return np.array([phi1, theta, phi2])


def euler_to_rotation_matrix(phi1, theta, phi2):
    """
    Convert Euler angles to rotation matrix (Bunge convention: Z-X-Z).

    Parameters:
    -----------
    phi1, theta, phi2 : float
        Euler angles in radians

    Returns:
    --------
    numpy.ndarray
        3x3 rotation matrix
    """
    cos_phi1, sin_phi1 = np.cos(phi1), np.sin(phi1)
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    cos_phi2, sin_phi2 = np.cos(phi2), np.sin(phi2)

    R = np.array(
        [
            [
                cos_phi1 * cos_phi2 - sin_phi1 * sin_phi2 * cos_theta,
                sin_phi1 * cos_phi2 + cos_phi1 * sin_phi2 * cos_theta,
                sin_phi2 * sin_theta,
            ],
            [
                -cos_phi1 * sin_phi2 - sin_phi1 * cos_phi2 * cos_theta,
                -sin_phi1 * sin_phi2 + cos_phi1 * cos_phi2 * cos_theta,
                cos_phi2 * sin_theta,
            ],
            [sin_phi1 * sin_theta, -cos_phi1 * sin_theta, cos_theta],
        ]
    )

    return R
