"""
This module contains the core scientific computing functions for the VPSC simulation.
"""

import numpy as np
from tensor import voigt, lu_inverse
from hardening import update_hardening
from orientation import update_orientation
from shape import update_shape
from twinning import update_twinning
from schmid import update_schmid_tensors


def chg_basis(*args, **kwargs):
    """
    Change tensor basis (stub for CHG_BASIS FORTRAN routine).
    TODO: Implement full logic as in FORTRAN.
    """
    raise NotImplementedError("chg_basis not yet implemented")


def calculate_elastic_moduli(config, step):
    """
    Calculates upper-bound (UB) or self-consistent (SC) elastic moduli.
    """
    print(f"Calculating elastic moduli (step {step})...")
    # Placeholder
    return config


def run_elastic_calculation(config):
    """
    Performs a thermo-elastic self-consistent calculation.
    """
    print("Running thermo-elastic calculation...")
    # Placeholder: returns dummy stress and strain rate tensors
    sbar = np.zeros(6)
    dbar = np.zeros(6)
    return sbar, dbar


def run_vpfc(config, step):
    """
    Performs a visco-plastic full-constraint (Taylor) calculation.
    """
    print(f"Running VPFC (Taylor) calculation for step {step}...")

    # Get applied boundary conditions
    applied_udot = config.get("applied_udot", np.zeros((3, 3)))
    iudot = config.get("iudot", np.ones((3, 3), dtype=int))
    eqincr = config.get("eqincr", 0.001)

    # Apply incremental strain based on boundary conditions
    # Convert 3x3 strain rate to 6-component Voigt notation
    # Voigt: [11, 22, 33, 23, 13, 12]
    dbar = np.zeros(6)
    dbar[0] = applied_udot[0, 0] * eqincr  # e11
    dbar[1] = applied_udot[1, 1] * eqincr  # e22
    dbar[2] = applied_udot[2, 2] * eqincr  # e33
    dbar[3] = applied_udot[1, 2] * eqincr  # e23 (gamma23/2)
    dbar[4] = applied_udot[0, 2] * eqincr  # e13 (gamma13/2)
    dbar[5] = applied_udot[0, 1] * eqincr  # e12 (gamma12/2)

    # For Taylor model, all grains have same strain rate
    # Calculate stress based on material properties
    grains = []
    sbar = np.zeros(6)
    total_weight = 0.0

    for iph, phase in enumerate(config.get("phases", [])):
        ngrains = phase.get("ngrains", 1)
        weight = 1.0 / ngrains if ngrains > 0 else 1.0

        # Get elastic constants for stress calculation
        elastic_matrix = phase.get("elastic_matrix", np.eye(6))

        for igr in range(ngrains):
            # For elastic-plastic, stress = C : strain_elastic
            # Simplified: assume elastic response for now
            grain_stress = elastic_matrix @ dbar
            sbar += weight * grain_stress
            total_weight += weight

    if total_weight > 0:
        sbar /= total_weight

    # Store current state for output
    config["current_strain"] = config.get("current_strain", np.zeros(6)) + dbar
    config["current_stress"] = sbar.copy()

    return sbar, dbar


def run_vpsc(config, step):
    """
    Performs a visco-plastic self-consistent (VPSC) calculation.
    """
    print(f"Running VPSC calculation for step {step}...")

    # Get applied boundary conditions
    applied_udot = config.get("applied_udot", np.zeros((3, 3)))
    iudot = config.get("iudot", np.ones((3, 3), dtype=int))
    eqincr = config.get("eqincr", 0.001)

    # Apply incremental strain based on boundary conditions
    # Convert 3x3 strain rate to 6-component Voigt notation
    dbar_applied = np.zeros(6)
    dbar_applied[0] = applied_udot[0, 0] * eqincr  # e11
    dbar_applied[1] = applied_udot[1, 1] * eqincr  # e22
    dbar_applied[2] = applied_udot[2, 2] * eqincr  # e33
    dbar_applied[3] = applied_udot[1, 2] * eqincr  # e23 (gamma23/2)
    dbar_applied[4] = applied_udot[0, 2] * eqincr  # e13 (gamma13/2)
    dbar_applied[5] = applied_udot[0, 1] * eqincr  # e12 (gamma12/2)

    max_outer = config.get("itmaxext", 50)
    tol = config.get("errs", 1e-5)
    nph = len(config.get("phases", []))

    # Initial guess for aggregate stress/strain
    sbar = np.zeros(6)
    dbar = dbar_applied.copy()  # Start with applied strain rate
    converged = False

    for outer in range(max_outer):
        sbar_old = sbar.copy()
        dbar_old = dbar.copy()

        # 1. Update effective medium properties (Eshelby, interaction tensor, etc.)
        #    (Placeholder: use identity tensors)
        C_eff = np.identity(6)

        # 2. For each grain, solve for local stress/strain given current medium
        sbar = np.zeros(6)
        dbar = dbar_applied.copy()  # Maintain applied boundary conditions
        total_weight = 0.0

        for iph, phase in enumerate(config.get("phases", [])):
            ngrains = phase.get("ngrains", 1)
            weight = 1.0 / ngrains if ngrains > 0 else 1.0

            # Get elastic constants for stress calculation
            elastic_matrix = phase.get("elastic_matrix", np.eye(6))

            for igr in range(ngrains):
                # Calculate stress from applied strain
                grain_stress = elastic_matrix @ dbar_applied
                sbar += weight * grain_stress
                total_weight += weight

        if total_weight > 0:
            sbar /= total_weight

        # 3. Check convergence (change in sbar/dbar < tol)
        delta = np.linalg.norm(sbar - sbar_old) + np.linalg.norm(dbar - dbar_old)
        if delta < tol:
            converged = True
            print(f"VPSC converged in {outer+1} iterations. Delta={delta:.2e}")
            break

    if not converged:
        print("Warning: VPSC did not converge within the maximum number of iterations.")

    # Store current state for output
    config["current_strain"] = config.get("current_strain", np.zeros(6)) + dbar
    config["current_stress"] = sbar.copy()

    return sbar, dbar
