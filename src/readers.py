"""
This module handles reading all input files for the VPSC-Python simulation.
"""

import numpy as np
import os


def parse_line(line):
    """Splits a line into a list of values, converting to float/int where possible."""
    parts = line.strip().split()
    converted_parts = []
    for part in parts:
        try:
            converted_parts.append(int(part))
        except ValueError:
            try:
                converted_parts.append(float(part))
            except ValueError:
                converted_parts.append(part)
    return converted_parts


def read_vpsc_input(file_path, log_file):
    """
    Reads the main VPSC input file (vpsc8.in) and populates a configuration dictionary.
    """
    print(f"Reading VPSC input from: {file_path}")
    log_file.write(f"**** INPUT FILE {os.path.basename(file_path)} FOR THIS RUN ****\n")

    with open(file_path, "r") as f:
        lines = f.readlines()

    for line in lines:
        log_file.write(line)
    log_file.write("\n")

    line_iter = iter(lines)
    config = {}

    # --- Initializations ---
    config["pi"] = np.pi
    config["xid3"] = np.identity(3)
    config["xid5"] = np.identity(5)
    config["xid6"] = np.identity(6)

    # --- Read Simulation Parameters ---
    config["iregime"] = parse_line(next(line_iter))[0]
    if config["iregime"] == -1:
        config["interaction"] = -1

    config["nph"] = parse_line(next(line_iter))[0]
    config["wph"] = parse_line(next(line_iter))

    # Normalize phase fractions
    wph_total = sum(config["wph"])
    if abs(wph_total - 1.0) > 1e-5:
        print("--> vol fraction of phases should add to 1 !!")
        print("--> will normalize phase fractions to 1 !!")
        config["wph"] = [w / wph_total for w in config["wph"]]

    # --- Phase-specific Data ---
    config["phases"] = []
    for iph in range(config["nph"]):
        phase_data = {}
        next(line_iter)  # Skip blank line/prosa

        line = parse_line(next(line_iter))
        phase_data["ishape"] = line[0]
        phase_data["ifrag"] = line[1]
        phase_data["crit_shp"] = line[2]

        phase_data["axisph"] = parse_line(next(line_iter))
        phase_data["eulerph"] = parse_line(next(line_iter))

        next(line_iter)  # Skip prosa
        phase_data["filetext"] = next(line_iter).strip()
        next(line_iter)  # Skip prosa
        phase_data["filecrys"] = next(line_iter).strip()
        next(line_iter)  # Skip prosa
        phase_data["fileaxes"] = next(line_iter).strip()
        next(line_iter)  # Skip prosa
        phase_data["idiff"] = parse_line(next(line_iter))[0]
        next(line_iter)  # Skip prosa
        phase_data["filediff"] = next(line_iter).strip()

        # Read crystal and grain data from their respective files
        phase_data.update(read_crystal_data(phase_data["filecrys"]))
        phase_data.update(read_texture_data(phase_data["filetext"]))

        config["phases"].append(phase_data)

    # --- Convergence and I/O Settings ---
    if config.get("interaction") != -1:
        next(line_iter)  # Skip prosa
        errs = parse_line(next(line_iter))
        config["errs"], config["errd"], config["errm"], config["errso"] = errs

        itmax = parse_line(next(line_iter))
        config["itmaxext"], config["itmaxint"], config["itmaxso"] = itmax

        irsvar = parse_line(next(line_iter))
        config["irvar"], config["jxrsini"], config["jxrsfin"], config["jxrstep"] = (
            irsvar
        )

        next(line_iter)  # Skip prosa
        config["irecover"] = parse_line(next(line_iter))[0]
        config["isave"] = parse_line(next(line_iter))[0]
        config["icubcom"] = parse_line(next(line_iter))[0]
        config["nwrite"] = parse_line(next(line_iter))[0]

        # --- Modeling Conditions ---
        next(line_iter)  # Skip prosa
        inter = parse_line(next(line_iter))
        config["interaction"], config["neffgrx"] = inter

        iupd = parse_line(next(line_iter))
        config["iupdori"], config["iupdshp"], config["iupdhar"] = iupd

        config["nneigh"] = parse_line(next(line_iter))[0]
        config["iflu"] = parse_line(next(line_iter))[0]

    print("Finished reading VPSC input.")
    return config


def read_crystal_data(filecrys_path):
    """
    Placeholder for reading crystal data file.
    """
    print(f"  Reading crystal data from: {filecrys_path}")
    # In a real implementation, this would open and parse the file.
    return {"crystal_props": "..."}


def read_texture_data(filetext_path):
    """
    Placeholder for reading texture data file.
    """
    print(f"  Reading texture data from: {filetext_path}")
    # In a real implementation, this would open and parse the file.
    return {"texture_props": "..."}


def read_postmortem(file_path, config):
    """
    Reads the postmortem file to recover the state from a previous run.
    """
    print(f"Reading postmortem data from: {file_path}")
    # Placeholder implementation
    return config


def read_cubcomp(file_path, config, step):
    """
    Reads ideal rolling components from CUBCOMP.IN.
    """
    print(f"Reading cubcomp data from: {file_path}")
    # Placeholder implementation
    return config


def read_load_conditions(file_path, config):
    """
    Reads load conditions from the specified process file.
    """
    print(f"Reading load conditions from: {file_path}")
    # Placeholder implementation
    return config


def read_var_vel_grad(config):
    """
    Reads variable velocity gradient history.
    """
    print("Reading variable velocity gradient history")
    # Placeholder implementation
    return config
