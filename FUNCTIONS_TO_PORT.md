# VPSC Python Port: Functions and Subroutines to Implement

This checklist covers the main routines from the FORTRAN VPSC codebase that need to be ported to Python for a faithful, full-featured conversion.

---

## Core VPSC Algorithm & Workflow
- `VPSC_INPUT` (input parsing, already partially ported)
- `VPFC` (Taylor model, skeleton present)
- `VPSC` (self-consistent model, skeleton present)
- `RUN_ELASTIC_CALCULATION` / `ELSC` (thermo-elastic calculation)
- `CALCULATE_ELASTIC_MODULI` (already stubbed)
- `LOAD_CONDITIONS` (process/load history parsing)
- `POSTMORTEM` (state recovery/output)

---

## Grain/Phase Physics
- `GRAIN_STRESS` (core grain stress solver, Newton-Raphson, etc.)
- `GRAIN_RATE_AND_MODULI` (grain strain rate and moduli update)
- `GRAIN_INFO` (grain initialization)
- `DATA_CRYSTAL` (crystal properties, slip/twin systems)
- `DATA_GRAIN` (grain orientations, weights)
- `INITIAL_STATE_GUESS` (initialization for iterative solvers)

---

## Tensor Algebra & Basis Changes
- `VOIGT` (Voigt notation conversions, ported)
- `CHG_BASIS` (basis change, ported)
- `LU_INVERSE`, `LU_DECOMP`, `LU_BACKSUBS`, `LU_EQSYSTEM` (matrix algebra, partially ported)
- `EULER` (Euler angle to rotation matrix)
- `VNORM`, `TNORM`, `VMISMATCH`, `TMISMATCH` (norms and mismatch metrics)
- `DET` (determinant)
- `RANDOM2` (random number generator, if needed)

---

## Eshelby & Effective Medium
- `ESHELBY_TENSOR` (Eshelby tensor calculation)
- `ESHELBY` (Eshelby integrals)
- `ESH_GAUSS_LEGENDRE`, `ESH_INV3_VOIGT`, `ESH_INV4_VOIGT`, `ESH_MULT_VOIGT` (integration and tensor algebra for Eshelby)

---

## Hardening & Constitutive Laws
- `UPDATE_HARDENING` (already stubbed)
- `UPDATE_CRSS_DD`, `UPDATE_CRSS_DD_REV`, `UPDATE_CRSS_MTS`, `UPDATE_CRSS_VOCE` (various hardening laws)
- `CHECK_VOCE` (Voce law check)
- `SCALE_GAMD0S` (shear rate scaling)

---

## Microstructure Evolution
- `UPDATE_ORIENTATION` (already stubbed)
- `UPDATE_SHAPE` (already stubbed)
- `UPDATE_TWINNING_PTR`, `UPDATE_TWINNING_VFT`, `TWIN_ORIENTATION` (twinning logic)
- `UPDATE_GROWTH_RATE`, `UPDATE_GROWTH_ZR`, `UPDATE_GROWTH_ZIRC2` (grain growth, if needed)

---

## Schmid Tensors & Slip/Twin Systems
- `UPDATE_SCHMID` (already stubbed)
- `CRYSTAL_SYMMETRY` (symmetry operations)
- `CUBCOMP` (ideal rolling components)
- `DIFF_PLANES` (diffraction planes, if needed)

---

## Polycrystal Yield Surface & Probing
- `PCYS`, `PCYS_IT` (yield surface probing)
- `STATE_NxN`, `STATE_6x6` (state equation solvers)

---

## Statistical Output & Analysis
- `STAT_GRAIN_SHAPE`, `STAT_SHEAR_ACTIVITY`, `STAT_STRESS_STRAIN`, `STAT_TWINNING` (statistical outputs)
- `WRITE_SHEAR_ACTIVITY`, `WRITE_STRESS_STRAIN`, `WRITE_TEXTURE`, `WRITE_TWIN_STATS` (output routines)

---

## Other Utilities
- `NEIGHBOURS` (neighbor assignment for grains)
- `N_EFFECTIVE` (effective n calculation)
- `NEWTON_RAPHSON` (nonlinear solver, used in grain stress)
- `LANKFORD` (Lankford coefficient calculation)
- `VAR_VEL_GRAD` (variable velocity gradient history)
- `TEXTURE_ROTATION` (texture update)
- `SO_*` (second-order routines, if needed for advanced models)

---

**Note:**
- Some routines are already partially or fully ported (e.g., `voigt`, `chg_basis`, some input/output).
- Many routines are currently placeholders or stubs in the Python codebase and need full scientific implementation and validation.
- Some routines (e.g., advanced growth, second-order, or diffraction routines) may be optional depending on your modeling needs.

---

**Recommended Next Steps:**
1. Prioritize porting and testing the core physics routines: `GRAIN_STRESS`, `GRAIN_RATE_AND_MODULI`, `ESHELBY_TENSOR`, `NEWTON_RAPHSON`, and the main VPSC/VPFC loops.
2. Then implement the supporting tensor algebra, hardening, and microstructure evolution routines.
3. Validate each module with unit tests and compare results to the FORTRAN code.
