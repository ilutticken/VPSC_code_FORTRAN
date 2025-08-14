def read_var_vel_grad(*args, **kwargs):
    """Stub for reading variable velocity gradient. Returns None."""
    return None


def read_load_conditions(proc_file, config):
    """Read process file containing loading conditions and boundary conditions."""
    import os
    import numpy as np

    # Handle relative paths - if file doesn't exist, try relative to current directory
    original_file = proc_file
    if not os.path.exists(proc_file):
        # Try looking in parent directory
        parent_file = os.path.join("..", proc_file)
        if os.path.exists(parent_file):
            proc_file = parent_file
        else:
            # Try extracting just the filename
            filename = os.path.basename(proc_file)
            if os.path.exists(filename):
                proc_file = filename
            elif os.path.exists(os.path.join("..", filename)):
                proc_file = os.path.join("..", filename)
            else:
                # Try going up to workspace root and look for file
                workspace_paths = [
                    os.path.join("..", "..", "..", filename),
                    os.path.join("..", "..", filename),
                    os.path.join(
                        "..", "..", proc_file.replace("\\", os.sep).replace("/", os.sep)
                    ),
                    os.path.join(
                        "..",
                        "..",
                        "..",
                        proc_file.replace("\\", os.sep).replace("/", os.sep),
                    ),
                ]
                for test_path in workspace_paths:
                    if os.path.exists(test_path):
                        proc_file = test_path
                        break

    if not os.path.exists(proc_file):
        print(
            f"Warning: Process file {original_file} not found, using default conditions"
        )
        return config

    try:
        with open(proc_file, "r") as f:
            lines = [line.strip() for line in f.readlines()]

        # Remove comments and empty lines, but keep numerical data lines
        clean_lines = []
        for line in lines:
            # Skip empty lines and lines starting with *
            if not line or line.startswith("*"):
                continue
            # Remove inline comments (everything after | character)
            if "|" in line:
                line = line.split("|")[0].strip()
            # Only keep lines that have numerical content
            if line and any(c.isdigit() or c in ".-+" for c in line):
                clean_lines.append(line)

        # print(f"Debug: Found {len(clean_lines)} clean lines in {proc_file}")
        # if len(clean_lines) > 0:
        #     print(f"Debug: First few lines: {clean_lines[:3]}")

        # First line: nsteps, ictrl, eqincr, tempini, tempfin
        if clean_lines:
            parts = clean_lines[0].split()
            if len(parts) >= 5:
                config["nsteps"] = int(parts[0])
                config["ictrl"] = int(parts[1])
                config["eqincr"] = float(parts[2])
                config["tempini"] = float(parts[3])
                config["tempfin"] = float(parts[4])

        # Parse boundary conditions
        if len(clean_lines) >= 8:
            # Velocity gradient flags (3x3 matrix)
            iudot = np.zeros((3, 3), dtype=int)
            for i in range(3):
                parts = clean_lines[1 + i].split()
                for j in range(3):
                    if j < len(parts):
                        iudot[i, j] = int(parts[j])

            # Velocity gradient values (3x3 matrix)
            udot = np.zeros((3, 3))
            for i in range(3):
                parts = clean_lines[4 + i].split()
                for j in range(3):
                    if j < len(parts):
                        udot[i, j] = float(parts[j])

            config["iudot"] = iudot
            config["udot"] = udot

            # Store applied strain rate for this step
            config["applied_udot"] = udot.copy()

            print(f"Loaded process file: {proc_file}")
            print(
                f"  Steps: {config['nsteps']}, Control: {config.get('ictrl', 'unknown')}"
            )
            print(f"  Applied velocity gradient:")
            print(f"    Flags:\n{iudot}")
            print(f"    Values:\n{udot}")

    except Exception as e:
        print(f"Warning: Error reading process file {proc_file}: {e}")

    return config


def read_cubcomp(file_path, config, step=0):
    """Stub for reading CUBCOMP.IN. Returns config unchanged."""
    return config


def read_postmortem(file_path, config):
    """
    Read binary POSTMORT.IN file to recover grain states from previous run.

    The FORTRAN POSTMORTEM subroutine handles binary I/O of grain states.
    This Python version mimics the structure for reading POSTMORT.IN files.

    Binary format (as per vpsc8sub.for):
    - isaveprevious (integer)
    - DAV(1:5) (real array - average deformation rate)
    - SAV(1:5) (real array - average stress)
    - DBAR0(1:5) (real array - reference deformation rate)
    - MBARTG(1:5,1:5) (real matrix - tangent modulus)
    - LBARTG(1:5,1:5) (real matrix - compliance)
    - SG(1:5,ngr1:ngr2) (real array - grain stresses)

    If isaveprevious != 1, also reads:
    - EPSTOTc(1:3,1:3), epsvm, epsacu (total strain, von Mises, accumulated)
    - axisph, fijph arrays (phase axes and deformation gradients)
    - For each phase: grain orientations, weights, hardening states, twin data
    """
    import struct
    import numpy as np

    if not os.path.exists(file_path):
        print(f"POSTMORT.IN file not found: {file_path}")
        return config

    try:
        with open(file_path, "rb") as f:
            # Read isaveprevious
            isaveprevious = struct.unpack("i", f.read(4))[0]
            config["isaveprevious"] = isaveprevious

            # Read average arrays (DAV, SAV, DBAR0)
            dav = struct.unpack("5f", f.read(20))
            sav = struct.unpack("5f", f.read(20))
            dbar0 = struct.unpack("5f", f.read(20))

            config["dav"] = np.array(dav)
            config["sav"] = np.array(sav)
            config["dbar0"] = np.array(dbar0)

            # Read tangent modulus and compliance matrices (5x5)
            mbartg = np.zeros((5, 5))
            lbartg = np.zeros((5, 5))
            for i in range(5):
                row = struct.unpack("5f", f.read(20))
                mbartg[i, :] = row
            for i in range(5):
                row = struct.unpack("5f", f.read(20))
                lbartg[i, :] = row

            config["mbartg"] = mbartg
            config["lbartg"] = lbartg

            # Read grain stresses SG(1:5, ngr1:ngr2)
            nph = config.get("nph", 1)
            ngr_total = config.get("ngr_total", 1000)  # Default grain count
            sg = np.zeros((5, ngr_total))
            for kkk in range(ngr_total):
                stress = struct.unpack("5f", f.read(20))
                sg[:, kkk] = stress
            config["sg"] = sg

            # If isaveprevious == 1, we're done (quick restart)
            if isaveprevious == 1:
                print(f"Read POSTMORT.IN for quick restart (isaveprevious=1)")
                return config

            # Otherwise, read full grain states
            # Total strain tensor, von Mises strain, accumulated strain
            epstotc = np.zeros((3, 3))
            for i in range(3):
                row = struct.unpack("3f", f.read(12))
                epstotc[i, :] = row
            epsvm, epsacu = struct.unpack("2f", f.read(8))

            config["epstotc"] = epstotc
            config["epsvm"] = epsvm
            config["epsacu"] = epsacu

            # Phase axes and deformation gradients
            # axisph(0:3, 1:3, 0:nph), fijph(1:3, 1:3, 0:nph)
            axisph = np.zeros((4, 3, nph + 1))
            fijph = np.zeros((3, 3, nph + 1))

            # Read for phase 0 (aggregate)
            for i in range(4):
                row = struct.unpack("3f", f.read(12))
                axisph[i, :, 0] = row
            for i in range(3):
                row = struct.unpack("3f", f.read(12))
                fijph[i, :, 0] = row

            config["axisph"] = axisph
            config["fijph"] = fijph

            # Read grain data for each phase
            grain_data = {}
            for iph in range(1, nph + 1):
                ngr1 = config.get(f"ngr_{iph-1}", 0) + 1
                ngr2 = config.get(f"ngr_{iph}", 1000)
                ngr_phase = ngr2 - ngr1 + 1

                # Phase axes and deformation gradient
                for i in range(4):
                    row = struct.unpack("3f", f.read(12))
                    axisph[i, :, iph] = row
                for i in range(3):
                    row = struct.unpack("3f", f.read(12))
                    fijph[i, :, iph] = row

                # Grain orientations ag(1:3, 1:3, ngr1:ngr2)
                ag = np.zeros((3, 3, ngr_phase))
                for kkk in range(ngr_phase):
                    for i in range(3):
                        row = struct.unpack("3f", f.read(12))
                        ag[i, :, kkk] = row

                # Grain weights and total shear
                wgt = struct.unpack(f"{ngr_phase}f", f.read(4 * ngr_phase))
                gtotgr = struct.unpack(f"{ngr_phase}f", f.read(4 * ngr_phase))

                grain_data[f"phase_{iph}"] = {
                    "ag": ag,
                    "wgt": np.array(wgt),
                    "gtotgr": np.array(gtotgr),
                }

                # Hardening data (depends on hardening law)
                ihardlaw = config.get("ihardlaw", 0)
                nsyst = config.get(f"nsyst_{iph}", 12)  # Default for FCC

                if ihardlaw == 0:  # Voce law
                    crss = np.zeros((nsyst, ngr_phase))
                    for kkk in range(ngr_phase):
                        row = struct.unpack(f"{nsyst}f", f.read(4 * nsyst))
                        crss[:, kkk] = row
                    grain_data[f"phase_{iph}"]["crss"] = crss
                elif ihardlaw == 1:  # MTS law
                    taue = np.zeros((nsyst, ngr_phase))
                    for kkk in range(ngr_phase):
                        row = struct.unpack(f"{nsyst}f", f.read(4 * nsyst))
                        taue[:, kkk] = row
                    grain_data[f"phase_{iph}"]["taue"] = taue

                # Twin data if present
                ntwmod = config.get(f"ntwmod_{iph}", 0)
                if ntwmod > 0:
                    pritw, sectw = struct.unpack("2f", f.read(8))
                    ntwevents = struct.unpack(f"{ngr_phase}i", f.read(4 * ngr_phase))
                    twfrph = struct.unpack(f"{ntwmod}f", f.read(4 * ntwmod))
                    eftwfr = struct.unpack(f"{ntwmod}f", f.read(4 * ntwmod))
                    ktwsmx = struct.unpack(f"{ngr_phase}i", f.read(4 * ngr_phase))

                    twfrsy = np.zeros((nsyst, ngr_phase))
                    for kkk in range(ngr_phase):
                        row = struct.unpack(f"{nsyst}f", f.read(4 * nsyst))
                        twfrsy[:, kkk] = row

                    grain_data[f"phase_{iph}"]["twin_data"] = {
                        "pritw": pritw,
                        "sectw": sectw,
                        "ntwevents": np.array(ntwevents),
                        "twfrph": np.array(twfrph),
                        "eftwfr": np.array(eftwfr),
                        "ktwsmx": np.array(ktwsmx),
                        "twfrsy": twfrsy,
                    }

            config["grain_data"] = grain_data
            print(f"Successfully read POSTMORT.IN with full grain states")

    except (FileNotFoundError, struct.error, IOError) as e:
        print(f"Error reading POSTMORT.IN: {e}")
        print("Continuing without recovery data")

    return config


def read_vpsc_input(input_path, log_handle=None):
    """
    Parse a VPSC input file and return a config dictionary.
    Handles the structure seen in ex03 examples.
    """
    config = {}
    with open(input_path, "r") as f:
        lines = [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("*")
        ]
    idx = 0
    # VP/ELAST regime
    config["regime"] = int(lines[idx].split()[0])
    idx += 1
    # Number of phases
    config["nph"] = int(lines[idx].split()[0])
    idx += 1
    # Relative volume fractions
    wph_tokens = lines[idx].split()
    config["wph"] = [float(wph_tokens[i]) for i in range(config["nph"])]
    idx += 1
    # Skip phase info header
    if "PHASE" in lines[idx].upper():
        idx += 1

    # Process each phase
    for phase in range(config["nph"]):
        # Grain shape, fragmentation, aspect ratio for this phase
        grain_shape_tokens = lines[idx].split()
        config[f"grain_shape_ph{phase+1}"] = [
            int(grain_shape_tokens[i]) for i in range(3)
        ]
        idx += 1
        # Initial ellipsoid ratios for this phase
        ellipsoid_ratios_tokens = lines[idx].split()
        config[f"ellipsoid_ratios_ph{phase+1}"] = [
            float(ellipsoid_ratios_tokens[i]) for i in range(3)
        ]
        idx += 1
        # Initial Euler angles for this phase
        ellipsoid_euler_tokens = lines[idx].split()
        config[f"ellipsoid_euler_ph{phase+1}"] = [
            float(ellipsoid_euler_tokens[i]) for i in range(3)
        ]
        idx += 1
        # Texture file for this phase
        config[f"filetext_ph{phase+1}"] = lines[idx]
        idx += 1
        # Crystal file for this phase
        config[f"filecrys_ph{phase+1}"] = lines[idx]
        idx += 1
        # Grain shape file for this phase
        config[f"fileaxes_ph{phase+1}"] = lines[idx]
        idx += 1
        # Diffraction file (may be 0 and dummy) for this phase
        diff_flag = lines[idx]
        idx += 1
        config[f"filediff_ph{phase+1}"] = lines[idx]
        idx += 1
        # Skip phase header for next phase
        if phase < config["nph"] - 1 and "PHASE" in lines[idx].upper():
            idx += 1
    # Precision settings
    precision_tokens = lines[idx].split()
    config["precision"] = [
        float(precision_tokens[i]) for i in range(min(4, len(precision_tokens)))
    ]
    idx += 1
    # Iteration limits
    itmax_tokens = lines[idx].split()
    config["itmax"] = [int(itmax_tokens[i]) for i in range(min(3, len(itmax_tokens)))]
    idx += 1
    # IRSVAR and related
    irsvar_tokens = lines[idx].split()
    config["irsvar"] = [
        int(irsvar_tokens[i]) for i in range(min(4, len(irsvar_tokens)))
    ]
    idx += 1
    # Input/output settings
    config["irecover"] = int(lines[idx].split()[0])
    idx += 1
    config["isave"] = int(lines[idx].split()[0])
    idx += 1
    config["icubcom"] = int(lines[idx].split()[0])
    idx += 1
    config["nwrite"] = int(lines[idx].split()[0])
    idx += 1
    # Modeling conditions
    interaction_tokens = lines[idx].split()
    config["interaction"] = int(interaction_tokens[0])
    # Check if second token is numeric (neffx) or comment
    if len(interaction_tokens) > 1:
        try:
            config["neffx"] = int(interaction_tokens[1])
        except ValueError:
            # Second token is a comment, use default neffx
            config["neffx"] = 10 if config["interaction"] == 3 else 1
    else:
        # Default neffx value when not specified
        config["neffx"] = 10 if config["interaction"] == 3 else 1
    idx += 1
    config["iupdate"] = [
        int(lines[idx].split()[i]) for i in range(min(3, len(lines[idx].split())))
    ]
    # Extract individual update flags for easier access
    iupdate = config["iupdate"]
    config["iupdori"] = iupdate[0] if len(iupdate) > 0 else 1  # Update orientation
    config["iupdshp"] = iupdate[1] if len(iupdate) > 1 else 1  # Update shape
    config["iupdhar"] = iupdate[2] if len(iupdate) > 2 else 1  # Update hardening
    idx += 1
    config["nneigh"] = int(lines[idx].split()[0])
    idx += 1
    config["iflu"] = int(lines[idx].split()[0])
    idx += 1
    # Number of processes
    config["nprocess"] = int(lines[idx].split()[0])
    idx += 1
    # Process blocks are handled in main.py
    if log_handle:
        print(f"Parsed VPSC input: {input_path}", file=log_handle)

    # Create phase data structures with crystal properties
    phases = []
    for phase_idx in range(config.get("nph", 1)):
        phase_data = {
            "phase_id": phase_idx + 1,
            "weight": (
                config["wph"][phase_idx]
                if phase_idx < len(config.get("wph", []))
                else 1.0
            ),
            "ngrains": 500,  # Default number of grains per phase
        }

        # Add grain shape parameters
        grain_shape_key = f"grain_shape_ph{phase_idx + 1}"
        if grain_shape_key in config:
            ishape, ifragm, crit_aspect = config[grain_shape_key]
            phase_data["ishape"] = ishape
            phase_data["ifragm"] = ifragm
            phase_data["crit_aspect_ratio"] = crit_aspect
        else:
            phase_data["ishape"] = 0  # No grain shape evolution
            phase_data["ifragm"] = 0  # No fragmentation
            phase_data["crit_aspect_ratio"] = 25  # Default

        # Add ellipsoid parameters
        ellipsoid_ratios_key = f"ellipsoid_ratios_ph{phase_idx + 1}"
        if ellipsoid_ratios_key in config:
            phase_data["ellipsoid_ratios"] = config[ellipsoid_ratios_key]
        else:
            phase_data["ellipsoid_ratios"] = [1.0, 1.0, 1.0]

        ellipsoid_euler_key = f"ellipsoid_euler_ph{phase_idx + 1}"
        if ellipsoid_euler_key in config:
            phase_data["ellipsoid_euler"] = config[ellipsoid_euler_key]
        else:
            phase_data["ellipsoid_euler"] = [0.0, 0.0, 0.0]

        # Load crystal properties if available
        filecrys_key = f"filecrys_ph{phase_idx + 1}"
        if filecrys_key in config:
            filecrys_path = config[filecrys_key]
            if not filecrys_path.lower() == "dummy":
                # Handle relative paths like we do for process files
                if not os.path.exists(filecrys_path):
                    # Try looking in parent directory
                    parent_file = os.path.join("..", filecrys_path)
                    if os.path.exists(parent_file):
                        filecrys_path = parent_file
                    else:
                        # Try extracting just the filename
                        filename = os.path.basename(filecrys_path)
                        if os.path.exists(filename):
                            filecrys_path = filename
                        elif os.path.exists(os.path.join("..", filename)):
                            filecrys_path = os.path.join("..", filename)

                crystal_data = read_crystal_data(filecrys_path)
                if "error" not in crystal_data:
                    # Use the already parsed elastic matrix from crystal data
                    if "elastic_matrix" in crystal_data:
                        phase_data["elastic_matrix"] = crystal_data["elastic_matrix"]
                        phase_data["crystal_symmetry"] = crystal_data.get(
                            "crysym", "cubic"
                        )

                        # Store growth parameters if present
                        phase_data["constitutive_law"] = crystal_data.get(
                            "constitutive_law", 0
                        )
                        if crystal_data.get("gamd0") is not None:
                            phase_data["gamd0"] = crystal_data["gamd0"]
                        if crystal_data.get("growth_rate_tensor") is not None:
                            phase_data["growth_rate_tensor"] = crystal_data[
                                "growth_rate_tensor"
                            ]

                        print(
                            f"  Loaded elastic constants for phase {phase_idx + 1}: {crystal_data.get('crysym', 'cubic')}"
                        )
                        if crystal_data.get("constitutive_law") == 30:
                            print(
                                f"  Irradiation growth model enabled (GAMD0={crystal_data.get('gamd0', 'N/A')})"
                            )
                    else:
                        print(
                            f"  Could not expand elastic constants for phase {phase_idx + 1}"
                        )
                        phase_data["elastic_matrix"] = (
                            np.eye(6) * 100e9
                        )  # Default stiffness
                else:
                    print(f"  Error loading crystal file for phase {phase_idx + 1}")
                    phase_data["elastic_matrix"] = (
                        np.eye(6) * 100e9
                    )  # Default stiffness
            else:
                phase_data["elastic_matrix"] = np.eye(6) * 100e9  # Default stiffness
        else:
            phase_data["elastic_matrix"] = np.eye(6) * 100e9  # Default stiffness

        phases.append(phase_data)

    config["phases"] = phases

    return config


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


def expand_elastic_constants(elastic, sym):
    """Expand elastic constants to 6x6 matrix for given symmetry."""
    import numpy as np

    C = np.zeros((6, 6))
    # CUBIC: 3 unique (C11, C12, C44)
    if sym.startswith("CUBIC") and len(elastic) >= 3:
        C11, C12, C44 = elastic[:3]
        C[:3, :3] = C12
        np.fill_diagonal(C[:3, :3], C11)
        C[3, 3] = C[4, 4] = C[5, 5] = C44
        return C
    # HEXAGONAL: 5 unique (C11, C12, C13, C33, C44)
    if sym.startswith("HEX") and len(elastic) >= 5:
        C11, C12, C13, C33, C44 = elastic[:5]
        C[0, 0] = C[1, 1] = C11
        C[0, 1] = C[1, 0] = C12
        C[0, 2] = C[2, 0] = C13
        C[1, 2] = C[2, 1] = C13
        C[2, 2] = C33
        C[3, 3] = C[4, 4] = C44
        C[5, 5] = 0.5 * (C11 - C12)
        return C
    # TETRAGONAL: 6 or 7 unique
    if sym.startswith("TETRA") and len(elastic) >= 6:
        # C11, C12, C13, C33, C44, C66
        C11, C12, C13, C33, C44, C66 = elastic[:6]
        C[0, 0] = C[1, 1] = C11
        C[0, 1] = C[1, 0] = C12
        C[0, 2] = C[2, 0] = C13
        C[1, 2] = C[2, 1] = C13
        C[2, 2] = C33
        C[3, 3] = C[4, 4] = C44
        C[5, 5] = C66
        return C
    # ORTHORHOMBIC: 9 unique
    if sym.startswith("ORTHO") and len(elastic) >= 9:
        # C11, C22, C33, C12, C13, C23, C44, C55, C66
        C11, C22, C33, C12, C13, C23, C44, C55, C66 = elastic[:9]
        C[0, 0] = C11
        C[1, 1] = C22
        C[2, 2] = C33
        C[0, 1] = C[1, 0] = C12
        C[0, 2] = C[2, 0] = C13
        C[1, 2] = C[2, 1] = C23
        C[3, 3] = C44
        C[4, 4] = C55
        C[5, 5] = C66
        return C
    # TRIGONAL: 6 unique
    if sym.startswith("TRIGON") and len(elastic) >= 6:
        # C11, C12, C13, C14, C33, C44
        C11, C12, C13, C14, C33, C44 = elastic[:6]
        C[0, 0] = C[1, 1] = C11
        C[0, 1] = C[1, 0] = C12
        C[0, 2] = C[2, 0] = C13
        C[1, 2] = C[2, 1] = C13
        C[2, 2] = C33
        C[0, 5] = C[5, 0] = C14
        C[1, 5] = C[5, 1] = -C14
        C[3, 3] = C[4, 4] = C44
        C[5, 5] = 0.5 * (C11 - C12)
        return C
    # MONOCLINIC: 13 unique


def read_crystal_data(filecrys_path):
    """
    Reads a .sx file (single crystal properties and slip/twin systems).
    """
    print(f"  Reading crystal data from: {filecrys_path}")
    data = {}
    if not os.path.exists(filecrys_path):
        print(f"  Crystal file not found: {filecrys_path}")
        return {"crystal_file": filecrys_path, "error": "not found"}
    with open(filecrys_path, "r") as f:
        # Keep both data lines and comment lines for constitutive law parsing
        all_lines = [line.strip() for line in f if line.strip()]
        # Regular lines (no comments) for elastic constants
        lines = [line.strip() for line in all_lines if not line.strip().startswith("*")]

    # Store raw lines for elastic constant parsing
    data["raw_lines"] = lines

    # Crystal symmetry and lattice
    # First line should be crystal symmetry (e.g., "cubic crysym")
    if len(lines) > 0 and "crysym" in lines[0]:
        data["crysym"] = (
            lines[0].split()[0].upper()
        )  # Extract "cubic" from "cubic crysym"
        data["material"] = "Unknown"  # Material name not in first line for this format
    else:
        data["material"] = lines[0] if len(lines) > 0 else "Unknown"
        data["crysym"] = lines[1].split()[0] if len(lines) > 1 else "cubic"

    data["lattice"] = (
        parse_line(lines[1]) if len(lines) > 1 else [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]
    )
    # Elastic constants (6x6, but may be incomplete or have comments)
    # Skip comment line "Elastic stiffness..." and start from data lines
    elastic = []
    start_line = 3  # Start after symmetry, lattice, and comment lines
    for i in range(start_line, min(start_line + 6, len(lines))):
        if i < len(lines):
            elastic.extend(
                [v for v in parse_line(lines[i]) if isinstance(v, (int, float, float))]
            )

    def expand_elastic_constants(elastic, sym):
        """Expand elastic constants to 6x6 matrix for given symmetry."""
        C = np.zeros((6, 6))
        # CUBIC: 3 unique (C11, C12, C44)
        if sym.startswith("CUBIC") and len(elastic) >= 3:
            C11, C12, C44 = elastic[:3]
            C[:3, :3] = C12
            np.fill_diagonal(C[:3, :3], C11)
            C[3, 3] = C[4, 4] = C[5, 5] = C44
            return C
        # HEXAGONAL: 5 unique (C11, C12, C13, C33, C44)
        if sym.startswith("HEX") and len(elastic) >= 5:
            C11, C12, C13, C33, C44 = elastic[:5]
            C[0, 0] = C[1, 1] = C11
            C[0, 1] = C[1, 0] = C12
            C[0, 2] = C[2, 0] = C13
            C[1, 2] = C[2, 1] = C13
            C[2, 2] = C33
            C[3, 3] = C[4, 4] = C44
            C[5, 5] = 0.5 * (C11 - C12)
            return C
        # TETRAGONAL: 6 or 7 unique
        if sym.startswith("TETRA") and len(elastic) >= 6:
            # C11, C12, C13, C33, C44, C66
            C11, C12, C13, C33, C44, C66 = elastic[:6]
            C[0, 0] = C[1, 1] = C11
            C[0, 1] = C[1, 0] = C12
            C[0, 2] = C[2, 0] = C13
            C[1, 2] = C[2, 1] = C13
            C[2, 2] = C33
            C[3, 3] = C[4, 4] = C44
            C[5, 5] = C66
            return C
        # ORTHORHOMBIC: 9 unique
        if sym.startswith("ORTHO") and len(elastic) >= 9:
            # C11, C22, C33, C12, C13, C23, C44, C55, C66
            C11, C22, C33, C12, C13, C23, C44, C55, C66 = elastic[:9]
            C[0, 0] = C11
            C[1, 1] = C22
            C[2, 2] = C33
            C[0, 1] = C[1, 0] = C12
            C[0, 2] = C[2, 0] = C13
            C[1, 2] = C[2, 1] = C23
            C[3, 3] = C44
            C[4, 4] = C55
            C[5, 5] = C66
            return C
        # TRIGONAL: 6 unique
        if sym.startswith("TRIGON") and len(elastic) >= 6:
            # C11, C12, C13, C14, C33, C44
            C11, C12, C13, C14, C33, C44 = elastic[:6]
            C[0, 0] = C[1, 1] = C11
            C[0, 1] = C[1, 0] = C12
            C[0, 2] = C[2, 0] = C13
            C[1, 2] = C[2, 1] = C13
            C[2, 2] = C33
            C[0, 5] = C[5, 0] = C14
            C[1, 5] = C[5, 1] = -C14
            C[3, 3] = C[4, 4] = C44
            C[5, 5] = 0.5 * (C11 - C12)
            return C
        # MONOCLINIC: 13 unique
        if sym.startswith("MONO") and len(elastic) >= 13:
            # C11, C22, C33, C44, C55, C66, C12, C13, C15, C23, C25, C35, C46
            C11, C22, C33, C44, C55, C66, C12, C13, C15, C23, C25, C35, C46 = elastic[
                :13
            ]
            C[0, 0] = C11
            C[1, 1] = C22
            C[2, 2] = C33
            C[3, 3] = C44
            C[4, 4] = C55
            C[5, 5] = C66
            C[0, 1] = C[1, 0] = C12
            C[0, 2] = C[2, 0] = C13
            C[0, 4] = C[4, 0] = C15
            C[1, 2] = C[2, 1] = C23
            C[1, 4] = C[4, 1] = C25
            C[2, 4] = C[4, 2] = C35
            C[3, 5] = C[5, 3] = C46
            return C
        # TRICLINIC: 21 unique
        if sym.startswith("TRICLIN") and len(elastic) >= 21:
            # Fill in row-wise
            vals = elastic[:21]
            idx = 0
            for row in range(6):
                for col in range(row, 6):
                    C[row, col] = C[col, row] = vals[idx]
                    idx += 1
            return C
        # If not recognized, return None
        return None

    # Thermal expansion (if present)
    for i in range(9, len(lines)):
        vals = [v for v in parse_line(lines[i]) if isinstance(v, (int, float, float))]
        if len(vals) >= 3:
            data["thermal_expansion"] = vals
            break

    # Create elastic matrix using crystal symmetry
    data["elastic_constants"] = elastic
    elastic_matrix = expand_elastic_constants(elastic, data["crysym"])
    if elastic_matrix is not None:
        data["elastic_matrix"] = elastic_matrix
        print(
            f"  Loaded {data['crysym']} crystal with elastic constants: {elastic[:3]}..."
        )
    else:
        print(
            f"  Warning: Could not parse elastic constants for {data['crysym']} symmetry"
        )
        data["elastic_matrix"] = np.eye(6)  # Fallback to identity

    # Parse slip/twin systems and constitutive law parameters
    data["slip_systems"] = []
    data["constitutive_law"] = 0  # Default Voce
    data["gamd0"] = None  # Irradiation creep parameter
    data["growth_rate_tensor"] = None  # Growth rate tensor

    # Search in all lines (including comments) for constitutive law
    slip_start = None
    constitutive_start = None
    for i, line in enumerate(all_lines):
        if "SLIP AND TWINNING MODES" in line:
            slip_start = i + 1
        elif "*Constitutive law" in line:
            # Found constitutive law comment - next line has the value
            if i + 1 < len(all_lines):
                try:
                    const_line = parse_line(all_lines[i + 1])
                    if const_line and len(const_line) >= 1:
                        data["constitutive_law"] = int(const_line[0])
                        print(f"  Found constitutive law: {data['constitutive_law']}")
                        constitutive_start = i + 1
                        break
                except:
                    pass

    # Parse growth parameters if constitutive law is 30
    if data["constitutive_law"] == 30:
        # Look for GAMD0 parameter in all lines
        for i, line in enumerate(all_lines):
            if "GAMD0" in line and i + 1 < len(all_lines):
                try:
                    gamd0_line = parse_line(all_lines[i + 1])
                    if gamd0_line and len(gamd0_line) >= 1:
                        data["gamd0"] = gamd0_line[0]
                        print(f"  Found GAMD0: {data['gamd0']}")
                        break
                except:
                    continue

        # Look for growth rate tensor in all lines
        for i, line in enumerate(all_lines):
            if "Growth rate tensor" in line:
                try:
                    # Read 3x3 tensor (next 3 lines)
                    tensor = []
                    for k in range(1, 4):
                        if i + k < len(all_lines):
                            row = parse_line(all_lines[i + k])
                            if row and len(row) >= 3:
                                tensor.append(row[:3])
                    if len(tensor) == 3:
                        data["growth_rate_tensor"] = np.array(tensor)
                        print(
                            f"  Found growth rate tensor: {data['growth_rate_tensor'][0,0]}, {data['growth_rate_tensor'][1,1]}, {data['growth_rate_tensor'][2,2]}"
                        )
                    break
                except:
                    continue

    # Basic slip system parsing (simplified)
    if slip_start:
        # This would need full implementation for complete .sx parsing
        # For now, we'll just note that slip systems are present
        data["has_slip_systems"] = True

    return data
