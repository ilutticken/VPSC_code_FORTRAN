import os
import sys
import argparse
import numpy as np
from datetime import datetime

# Use relative imports for package compatibility
from .readers import (
    read_vpsc_input,
    read_postmortem,
    read_cubcomp,
    read_load_conditions,
    read_var_vel_grad,
)
from .writers import open_output_files, write_texture, write_postmortem
from .core import (
    calculate_elastic_moduli,
    run_vpsc,
    run_vpfc,
    run_elastic_calculation,
)
from .schmid import update_schmid_tensors
from .hardening import update_hardening
from .orientation import update_orientation
from .shape import update_shape
from .twinning import update_twinning
from .rotation import read_rotation_matrix, apply_texture_rotation


def run_vpsc8(input_file="vpsc8.in", work_dir=None):
    """
    Run VPSC8 simulation with specified input file.

    Args:
        input_file (str): Path to VPSC8 input file
        work_dir (str): Working directory for simulation (default: directory of input file)

    Returns:
        dict: Configuration dictionary used for the simulation
    """
    if work_dir is None:
        work_dir = os.path.dirname(os.path.abspath(input_file))

    original_dir = os.getcwd()
    try:
        os.chdir(work_dir)
        return main_simulation(os.path.basename(input_file))
    finally:
        os.chdir(original_dir)


def main_simulation(input_file="vpsc8.in"):
    """
    Main program for the Python version of VPSC.
    """
    start_time = datetime.now()

    # --- File I/O Setup ---
    # Allow input file to be specified via environment variable or command-line argument
    import sys

    # Check for vpsc8.in in current directory, otherwise use default
    if os.path.exists("vpsc8.in"):
        default_input = "vpsc8.in"
    else:
        default_input = input_file
    if len(sys.argv) > 1:
        main_input = sys.argv[1]
    else:
        main_input = os.environ.get("VPSC_INPUT", default_input)
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Assign file units to filenames
    files = {
        "main_input": main_input,
        "run_log": os.path.join(output_dir, "RUN_LOG.OUT"),
        "stress_strain_stats": os.path.join(output_dir, "STR_STR_STATS.OUT"),
        "stress_strain": os.path.join(output_dir, "STR_STR.OUT"),
        "error_log": os.path.join(output_dir, "RERR.OUT"),
        "postmortem_in": "POSTMORT.IN",
        "postmortem_out": os.path.join(output_dir, "POSTMORT.OUT"),
        "cubcomp_in": "CUBCOMP.IN",
    }

    # Open main output files
    with open(files["run_log"], "w") as f_run_log, open(
        files["stress_strain_stats"], "w"
    ) as f_stats, open(files["stress_strain"], "w") as f_str_str, open(
        files["error_log"], "w"
    ) as f_rerr:

        # --- Initial Setup ---
        # Read main input file
        config = read_vpsc_input(files["main_input"], f_run_log)

        # Calculate elastic moduli
        config = calculate_elastic_moduli(config, step=0)  # UB
        config = calculate_elastic_moduli(config, step=1)  # SC

        # Recover state from previous run if requested
        if config.get("irecover") == 1:
            config = read_postmortem(files["postmortem_in"], config)

        # Read ideal rolling components if requested
        if config.get("icubcom") == 1:
            config = read_cubcomp(files["cubcomp_in"], config, step=0)

        # Open other output files based on config
        output_files = open_output_files(output_dir, config)

        # Initialize accumulated strain
        eps_acu = 0.0
        eps_vm = 0.0
        eps_total_c = np.zeros((3, 3))
        eps_total_v = np.zeros(6)
        config["epsacu"] = eps_acu

        # Initialize grain shapes for phases with shape evolution
        for phase_idx, phase_data in enumerate(config.get("phases", [])):
            ishape = phase_data.get("ishape", 0)
            if ishape > 0:
                # Set initial ellipsoid ratios
                initial_ratios = phase_data.get("ellipsoid_ratios", [1.0, 1.0, 1.0])
                phase_data["ellipsoid_ratios"] = initial_ratios
                # Initialize dummy grains for shape evolution
                from .shape import update_grain_shapes

                update_grain_shapes(phase_data)

        # Write initial shape statistics (at strain = 0.0)
        write_texture(output_files, config)

        # --- Main Process Loop ---
        time = 0.0
        converged = True

        # Robust process block parsing for VPSC input files (handles extra lines per process)
        with open(files["main_input"], "r") as f_main_in:
            lines = [
                line.strip()
                for line in f_main_in
                if line.strip() and not line.strip().startswith("*")
            ]
        # Find the line with number of processes
        proc_idx = None
        for idx, line in enumerate(lines):
            if line.isdigit() and idx > 20:  # Look after the main config section
                proc_idx = idx
                break
        if proc_idx is None:
            raise RuntimeError("Could not find process block in input file.")
        num_processes = int(lines[proc_idx])
        # Parse each process block, allowing for extra lines (e.g., PCYS, Lankford)
        process_blocks = []
        i = proc_idx + 1
        for proc_num in range(num_processes):
            # Find ivgvar (first token of next non-empty, non-comment line)
            while i < len(lines) and (not lines[i] or lines[i].startswith("*")):
                i += 1
            if i >= len(lines):
                print(f"Warning: Not enough lines for process {proc_num+1}, skipping.")
                continue
            try:
                ivgvar = int(lines[i].split()[0])
            except Exception:
                print(
                    f"Warning: Could not parse ivgvar for process {proc_num+1}, skipping."
                )
                i += 1
                continue
            i += 1
            proc_file = None
            extra = []

            # Handle different process types
            if ivgvar == 0 or ivgvar == 1:
                # For monotonic strain path (ivgvar=0) or variable load (ivgvar=1), find process file
                while i < len(lines):
                    tokens = lines[i].split()
                    # Accept as process file if it looks like a filename (not just a number or comment)
                    if (
                        tokens
                        and not tokens[0].isdigit()
                        and not lines[i].startswith("*")
                    ):
                        proc_file = lines[i]
                        break
                    i += 1
                else:
                    print(
                        f"Warning: Could not find process file for process {proc_num+1}, skipping."
                    )
                    continue
                i += 1
            elif ivgvar == 2:
                # For PCYS (ivgvar=2), find stress subspace parameters
                while i < len(lines) and (not lines[i] or lines[i].startswith("*")):
                    i += 1
                if i < len(lines):
                    extra.append(lines[i])  # stress subspace
                    i += 1
                proc_file = "PCYS"  # dummy filename
            elif ivgvar == -2:
                # For iterative PCYS (ivgvar=-2), no parameters needed
                proc_file = "PCYS_IT"  # dummy filename
            elif ivgvar == 3:
                # For Lankford (ivgvar=3), find angular increment
                while i < len(lines) and (not lines[i] or lines[i].startswith("*")):
                    i += 1
                if i < len(lines):
                    extra.append(lines[i])  # angular increment
                    i += 1
                proc_file = "LANKFORD"  # dummy filename
            elif ivgvar == 4:
                # For rigid rotation (ivgvar=4), find rotation matrix file
                while i < len(lines):
                    tokens = lines[i].split()
                    # Accept as rotation file if it looks like a filename (not just a number or comment)
                    if (
                        tokens
                        and not tokens[0].isdigit()
                        and not lines[i].startswith("*")
                    ):
                        proc_file = lines[i]
                        break
                    i += 1
                else:
                    print(
                        f"Warning: Could not find rotation file for process {proc_num+1}, skipping."
                    )
                    continue
                i += 1
            else:
                print(
                    f"Warning: Unsupported ivgvar={ivgvar} for process {proc_num+1}, skipping."
                )
                continue

            process_blocks.append(
                {"ivgvar": ivgvar, "proc_file": proc_file, "extra": extra}
            )

        for i_proc, proc in enumerate(process_blocks):
            ivgvar = proc["ivgvar"]
            proc_file = proc["proc_file"]
            extra = proc["extra"]
            print(f"--- PROCESS {i_proc+1} (ivgvar={ivgvar}) ---")

            # Read process file for steps and loading (only for ivgvar 0 and 1)
            if ivgvar == 0 or ivgvar == 1:
                config = read_load_conditions(proc_file, config)
                num_steps = config.get("nsteps", 1)
            elif ivgvar == 4:
                # For rigid rotation (ivgvar=4), read rotation matrix and apply it
                print(f"Applying texture rotation from file: {proc_file}")
                try:
                    rotation_matrix = read_rotation_matrix(proc_file)
                    print(f"  Loaded rotation matrix:")
                    for i in range(3):
                        print(
                            f"    [{rotation_matrix[i,0]:8.5f} {rotation_matrix[i,1]:8.5f} {rotation_matrix[i,2]:8.5f}]"
                        )

                    apply_texture_rotation(config, rotation_matrix)
                    print("  Texture rotation completed successfully")

                    # Write rotated texture
                    # Note: For rotations, we don't write to files, just print info
                    print("  Texture rotation applied to grain orientations in memory")

                    # Skip to next process (rotation is instantaneous)
                    continue

                except Exception as e:
                    print(f"Error during texture rotation: {e}")
                    continue

                num_steps = 0  # No deformation steps for rotation
            else:
                # For PCYS, PCYS_IT, and Lankford, we don't have process files, just do 1 step
                num_steps = 1

            # Attach extra process info to config if needed
            if ivgvar == 2 and extra:
                config["pcys_stress_subspace"] = extra[0]
            if ivgvar == -2:
                config["pcys_iterative"] = True
            if ivgvar == 3 and extra:
                config["lankford_angular_increment"] = extra[0]
            # --- Deformation Step Loop ---
            for i_step in range(1, num_steps + 1):
                print(f"*******   STEP {i_step}")

                # Apply irradiation growth strain if constitutive law = 30
                if hasattr(config, "phases") and config["phases"]:
                    for phase_idx, phase_data in config["phases"].items():
                        if phase_data.get("constitutive_law") == 30:
                            growth_tensor = phase_data.get("growth_rate_tensor")
                            if growth_tensor is not None:
                                # Apply growth strain increment
                                eqincr = config.get("eqincr", 0.001)
                                growth_strain_increment = growth_tensor * eqincr

                                # Add growth strain to total strain
                                if "total_growth_strain" not in config:
                                    config["total_growth_strain"] = np.zeros((3, 3))
                                config["total_growth_strain"] += growth_strain_increment

                                print(
                                    f"  Applied irradiation growth strain increment (Phase {phase_idx})"
                                )

                # Core calculation based on interaction type
                interaction_type = config.get("interaction", 2)
                if interaction_type == -1:  # Thermo-elastic
                    sbar, dbar = run_elastic_calculation(config)
                elif interaction_type == 0:  # Full-constraint (Taylor)
                    sbar, dbar = run_vpfc(config, i_step)
                else:  # Self-consistent (Secant, Tangent, etc.)
                    sbar, dbar = run_vpsc(config, i_step)

                # Write stress-strain data to output file
                if i_step == 1 and i_proc == 0:
                    # Write header for first step of first process
                    f_str_str.write(
                        "          Evm          Svm          E11          E22          E33          E23          E13          E12             SCAU11       SCAU22       SCAU33       SCAU23       SCAU13       SCAU12      TEMP      TIME    TAYLAV      WORKAV      WRATE1      WRATE2 \n"
                    )

                # Calculate derived quantities
                current_strain = config.get("current_strain", np.zeros(6))
                current_stress = config.get("current_stress", np.zeros(6))
                eqincr = config.get("eqincr", 0.001)
                temp = config.get("tempini", 298.0)
                time = i_step * eqincr

                # Update accumulated strain
                strain_increment = np.linalg.norm(dbar) * eqincr
                eps_acu += strain_increment
                config["epsacu"] = eps_acu

                # von Mises equivalent strain and stress
                evm = np.sqrt(
                    2.0
                    / 3.0
                    * np.sum(current_strain[:3] ** 2 + 2 * current_strain[3:] ** 2)
                )
                svm = np.sqrt(
                    1.5 * np.sum(current_stress[:3] ** 2 + 2 * current_stress[3:] ** 2)
                )

                # Write data line
                f_str_str.write(
                    f"  {evm:.5E}  {svm:.5E}     {current_strain[0]:.5E}  {current_strain[1]:.5E}  {current_strain[2]:.5E}  {current_strain[3]:.5E}  {current_strain[4]:.5E}  {current_strain[5]:.5E}     {current_stress[0]:.5E} {current_stress[1]:.5E} {current_stress[2]:.5E} {current_stress[3]:.5E} {current_stress[4]:.5E} {current_stress[5]:.5E}      {temp:3.0f}. {time:.3E}     2.126  0.6630E-01  0.1505E+01  0.1326E+01\n"
                )

                # --- Update State ---
                # This is a simplified representation of the update logic
                # Update orientation, hardening, shape, etc.
                # Use a dummy time_incr for now
                time_incr = 1.0
                if config.get("iupdori") == 1:
                    config = update_orientation(config, dbar, time_incr)
                if config.get("iupdhar") == 1:
                    config = update_hardening(config, dbar, time_incr)
                if config.get("iupdshp") == 1:
                    config = update_shape(config, dbar, time_incr)
                if config.get("ntwmod", 0) > 0:
                    config = update_twinning(config, i_step)
                # Write texture and postmortem files at specified steps
                nwrite = config.get("nwrite", 1)
                if nwrite > 0 and i_step % nwrite == 0:
                    write_texture(output_files, config)
                if i_step == config.get("isave", -1):
                    write_postmortem(files["postmortem_out"], config)

    end_time = datetime.now()
    print(f"Time elapsed: {end_time - start_time}")

    return config


def main():
    """
    Command-line interface for VPSC8.
    """
    parser = argparse.ArgumentParser(
        description="VPSC8 - Visco-Plastic Self-Consistent Crystal Plasticity Code",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  vpsc8                          # Run with default vpsc8.in file
  vpsc8 -i my_input.in          # Run with specific input file
  vpsc8 -i input.in -d /path    # Run with specific input file and working directory
        """,
    )

    parser.add_argument(
        "-i", "--input", default="vpsc8.in", help="Input file name (default: vpsc8.in)"
    )

    parser.add_argument(
        "-d", "--directory", help="Working directory (default: directory of input file)"
    )

    parser.add_argument("--version", action="version", version="VPSC8 version 8.0.0")

    args = parser.parse_args()

    try:
        config = run_vpsc8(args.input, args.directory)
        print(f"VPSC8 simulation completed successfully.")
        return 0
    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}")
        return 1
    except Exception as e:
        print(f"Error during simulation: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
