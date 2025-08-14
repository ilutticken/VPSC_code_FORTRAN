import numpy as np


def update_shape(config, dbar, time_incr):
    """
    Updates the shape of each grain based on ishape setting.

    Parameters:
    config: Configuration dictionary containing phase information
    dbar: Deformation rate tensor (strain rate)
    time_incr: Time increment for this step

    Returns:
    config: Updated configuration with new grain shapes
    """

    # Check if any phases have grain shape evolution enabled
    has_shape_evolution = False
    if "phases" in config:
        for phase_data in config["phases"]:
            if phase_data.get("ishape", 0) > 0:
                has_shape_evolution = True
                break

    if not has_shape_evolution:
        return config

    print("Updating grain shape...")

    # Process each phase
    for phase_idx, phase_data in enumerate(config.get("phases", [])):
        ishape = phase_data.get("ishape", 0)

        if ishape == 0:
            # No grain shape evolution
            continue

        elif ishape == 1:
            # Shape evolution via principal plastic velocity gradients
            # Use the plastic deformation rate to update ellipsoid axes
            update_shape_plastic_rate(phase_data, dbar, time_incr)
            # Also update individual grains with the same ellipsoid ratios
            update_grain_shapes(phase_data)

        elif ishape == 2:
            # Shape evolution via deformation gradient
            # Use total deformation gradient to update ellipsoid
            update_shape_deformation_gradient(phase_data, dbar, time_incr)
            # Also update individual grains with the same ellipsoid ratios
            update_grain_shapes(phase_data)

        elif ishape == 3:
            # Shape prescribed from file (not implemented in this version)
            print(f"  Phase {phase_idx+1}: Shape from file (not implemented)")

        elif ishape == 4:
            # Spherical grains (no evolution)
            phase_data["ellipsoid_ratios"] = [1.0, 1.0, 1.0]

    return config


def update_grain_shapes(phase_data):
    """
    Update individual grain ellipsoid ratios from phase-level data.
    """
    ellipsoid_ratios = phase_data.get("ellipsoid_ratios", [1.0, 1.0, 1.0])

    # If grains don't exist yet, create some dummy grains for shape evolution
    if "grains" not in phase_data or not phase_data["grains"]:
        import numpy as np

        # Create dummy grains for testing
        ngrains = 500  # Standard number
        phase_data["grains"] = []
        for i in range(ngrains):
            # Add small random variations to ellipsoid ratios for realistic statistics
            variations = np.random.normal(0, 0.01, 3)  # Small random noise
            grain_ratios = [
                max(0.1, ellipsoid_ratios[j] + variations[j]) for j in range(3)
            ]
            phase_data["grains"].append(
                {"ellipsoid_ratios": grain_ratios, "grain_id": i}
            )
    else:
        # Update all grains in this phase
        for grain in phase_data.get("grains", []):
            grain["ellipsoid_ratios"] = ellipsoid_ratios.copy()


def update_shape_plastic_rate(phase_data, dbar, time_incr):
    """
    Update grain shape using plastic velocity gradients (ishape=1).
    """
    # Get current ellipsoid ratios
    current_ratios = phase_data.get("ellipsoid_ratios", [1.0, 1.0, 1.0])

    # Convert strain rate tensor to principal values
    # dbar is assumed to be in Voigt notation [d11, d22, d33, d23, d13, d12]
    if len(dbar) == 6:
        # Convert from Voigt to matrix form
        dbar_matrix = np.array(
            [
                [dbar[0], dbar[5], dbar[4]],
                [dbar[5], dbar[1], dbar[3]],
                [dbar[4], dbar[3], dbar[2]],
            ]
        )
    else:
        dbar_matrix = dbar

    # Find principal strain rates
    eigenvals, eigenvecs = np.linalg.eigh(dbar_matrix)

    # Sort by magnitude (largest first)
    idx = np.argsort(np.abs(eigenvals))[::-1]
    principal_rates = eigenvals[idx]

    # Update ellipsoid ratios based on principal strain rates
    # Integrate: d(ln(a))/dt = d_principal
    dt = time_incr
    for i in range(3):
        if i < len(current_ratios):
            # Semi-logarithmic integration
            current_ratios[i] *= np.exp(principal_rates[i] * dt)

    # Normalize to maintain volume (for plastic deformation)
    volume_factor = np.prod(current_ratios) ** (1 / 3)
    current_ratios = [r / volume_factor for r in current_ratios]

    # Apply aspect ratio limits
    crit_aspect = phase_data.get("crit_aspect_ratio", 25.0)
    max_ratio = max(current_ratios)
    min_ratio = min(current_ratios)
    if max_ratio / min_ratio > crit_aspect:
        # Rescale to maintain critical aspect ratio
        scale_factor = (crit_aspect * min_ratio / max_ratio) ** 0.5
        current_ratios = [
            r * scale_factor if r == max_ratio else r / scale_factor
            for r in current_ratios
        ]

    phase_data["ellipsoid_ratios"] = current_ratios


def update_shape_deformation_gradient(phase_data, dbar, time_incr):
    """
    Update grain shape using deformation gradient (ishape=2).
    """
    # Get current ellipsoid ratios
    current_ratios = phase_data.get("ellipsoid_ratios", [1.0, 1.0, 1.0])

    # For ishape=2, we use the velocity gradient directly
    # Convert strain rate to velocity gradient (assuming small strains)
    if len(dbar) == 6:
        # Convert from Voigt to matrix form
        dbar_matrix = np.array(
            [
                [dbar[0], dbar[5], dbar[4]],
                [dbar[5], dbar[1], dbar[3]],
                [dbar[4], dbar[3], dbar[2]],
            ]
        )
    else:
        dbar_matrix = dbar

    # Integrate deformation gradient: F = exp(L*dt)
    dt = time_incr
    L = dbar_matrix  # Velocity gradient (symmetric part)

    # For small time steps, use first-order approximation
    F_increment = np.eye(3) + L * dt

    # Get current deformation gradient (start with identity if not stored)
    if "deformation_gradient" not in phase_data:
        phase_data["deformation_gradient"] = np.eye(3)

    F_total = F_increment @ phase_data["deformation_gradient"]
    phase_data["deformation_gradient"] = F_total

    # Extract ellipsoid evolution from deformation gradient
    # Use singular value decomposition: F = U * S * V^T
    U, S, Vt = np.linalg.svd(F_total)

    # The singular values give the principal stretches
    principal_stretches = S

    # Update ellipsoid ratios (normalize to maintain relative shape)
    mean_stretch = np.mean(principal_stretches)
    current_ratios = (principal_stretches / mean_stretch).tolist()

    # Apply aspect ratio limits
    crit_aspect = phase_data.get("crit_aspect_ratio", 25.0)
    max_ratio = max(current_ratios)
    min_ratio = min(current_ratios)
    if max_ratio / min_ratio > crit_aspect:
        # Rescale to maintain critical aspect ratio
        scale_factor = np.sqrt(crit_aspect * min_ratio / max_ratio)
        current_ratios = [
            r * scale_factor if r == max_ratio else r / scale_factor
            for r in current_ratios
        ]

    phase_data["ellipsoid_ratios"] = current_ratios
