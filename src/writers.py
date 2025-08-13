"""
This module handles writing all output files for the VPSC-Python simulation.
"""

import os


def open_output_files(output_dir, config):
    """
    Opens various output files based on the simulation configuration.
    """
    print("Opening output files...")
    files = {}
    nfiles = config.get("nph", 1)
    flabel = "PH"  # or 'EL' depending on nelem

    for i in range(1, nfiles + 1):
        files[f"tex_{i}"] = open(os.path.join(output_dir, f"TEX_{flabel}{i}.OUT"), "w")
        # Get ishape from the phase dictionary if available
        phase_ishape = 0
        if "phases" in config and len(config["phases"]) >= i:
            phase_ishape = config["phases"][i - 1].get("ishape", 0)
        if phase_ishape != 0:
            files[f"mor_{i}"] = open(
                os.path.join(output_dir, f"MOR_{flabel}{i}.OUT"), "w"
            )
            files[f"stat_axes_{i}"] = open(
                os.path.join(output_dir, f"STAT_AXES_{flabel}{i}.OUT"), "w"
            )
        if config.get("icubcom") == 1:
            files[f"cubcomp_{i}"] = open(
                os.path.join(output_dir, f"CUBCOMP{i}.OUT"), "w"
            )
        files[f"act_{i}"] = open(os.path.join(output_dir, f"ACT_{flabel}{i}.OUT"), "w")

    return files


def write_texture(output_files, config):
    """
    Writes the texture files.
    """
    print("Writing texture files...")
    # Placeholder implementation
    pass


def write_postmortem(file_path, config):
    """
    Write binary POSTMORT.OUT file to save grain states for future recovery.

    This Python version mimics the FORTRAN POSTMORTEM subroutine for writing
    grain states to a binary file that can be read by read_postmortem().

    Binary format (as per vpsc8sub.for):
    - isave (integer)
    - DAV(1:5) (real array - average deformation rate)
    - SAV(1:5) (real array - average stress)
    - DBAR0(1:5) (real array - reference deformation rate)
    - MBARTG(1:5,1:5) (real matrix - tangent modulus)
    - LBARTG(1:5,1:5) (real matrix - compliance)
    - SG(1:5,ngr1:ngr2) (real array - grain stresses)

    If isave != 1, also writes:
    - EPSTOTc(1:3,1:3), epsvm, epsacu (total strain, von Mises, accumulated)
    - axisph, fijph arrays (phase axes and deformation gradients)
    - For each phase: grain orientations, weights, hardening states, twin data
    """
    import struct
    import numpy as np

    print(f"Writing postmortem data to: {file_path}")

    try:
        with open(file_path, "wb") as f:
            # Write isave
            isave = config.get("isave", 1)
            f.write(struct.pack("i", isave))

            # Write average arrays (DAV, SAV, DBAR0)
            dav = config.get("dav", np.zeros(5))
            sav = config.get("sav", np.zeros(5))
            dbar0 = config.get("dbar0", np.zeros(5))

            f.write(struct.pack("5f", *dav))
            f.write(struct.pack("5f", *sav))
            f.write(struct.pack("5f", *dbar0))

            # Write tangent modulus and compliance matrices (5x5)
            mbartg = config.get("mbartg", np.eye(5))
            lbartg = config.get("lbartg", np.eye(5))

            for i in range(5):
                f.write(struct.pack("5f", *mbartg[i, :]))
            for i in range(5):
                f.write(struct.pack("5f", *lbartg[i, :]))

            # Write grain stresses SG(1:5, ngr1:ngr2)
            nph = config.get("nph", 1)
            ngr_total = config.get("ngr_total", 1000)
            sg = config.get("sg", np.zeros((5, ngr_total)))

            for kkk in range(ngr_total):
                f.write(struct.pack("5f", *sg[:, kkk]))

            # If isave == 1, we're done (quick save)
            if isave == 1:
                print(f"Wrote POSTMORT.OUT for quick restart (isave=1)")
                return

            # Otherwise, write full grain states
            # Total strain tensor, von Mises strain, accumulated strain
            epstotc = config.get("epstotc", np.zeros((3, 3)))
            epsvm = config.get("epsvm", 0.0)
            epsacu = config.get("epsacu", 0.0)

            for i in range(3):
                f.write(struct.pack("3f", *epstotc[i, :]))
            f.write(struct.pack("2f", epsvm, epsacu))

            # Phase axes and deformation gradients
            axisph = config.get("axisph", np.zeros((4, 3, nph + 1)))
            fijph = config.get("fijph", np.zeros((3, 3, nph + 1)))

            # Write for phase 0 (aggregate)
            for i in range(4):
                f.write(struct.pack("3f", *axisph[i, :, 0]))
            for i in range(3):
                f.write(struct.pack("3f", *fijph[i, :, 0]))

            # Write grain data for each phase
            grain_data = config.get("grain_data", {})
            for iph in range(1, nph + 1):
                ngr1 = config.get(f"ngr_{iph-1}", 0) + 1
                ngr2 = config.get(f"ngr_{iph}", 1000)
                ngr_phase = ngr2 - ngr1 + 1

                # Phase axes and deformation gradient
                for i in range(4):
                    f.write(struct.pack("3f", *axisph[i, :, iph]))
                for i in range(3):
                    f.write(struct.pack("3f", *fijph[i, :, iph]))

                # Get grain data for this phase
                phase_data = grain_data.get(f"phase_{iph}", {})

                # Grain orientations ag(1:3, 1:3, ngr1:ngr2)
                ag = phase_data.get("ag", np.zeros((3, 3, ngr_phase)))
                for kkk in range(ngr_phase):
                    for i in range(3):
                        f.write(struct.pack("3f", *ag[i, :, kkk]))

                # Grain weights and total shear
                wgt = phase_data.get("wgt", np.ones(ngr_phase))
                gtotgr = phase_data.get("gtotgr", np.zeros(ngr_phase))

                f.write(struct.pack(f"{ngr_phase}f", *wgt))
                f.write(struct.pack(f"{ngr_phase}f", *gtotgr))

                # Hardening data (depends on hardening law)
                ihardlaw = config.get("ihardlaw", 0)
                nsyst = config.get(f"nsyst_{iph}", 12)  # Default for FCC

                if ihardlaw == 0:  # Voce law
                    crss = phase_data.get("crss", np.ones((nsyst, ngr_phase)))
                    for kkk in range(ngr_phase):
                        f.write(struct.pack(f"{nsyst}f", *crss[:, kkk]))
                elif ihardlaw == 1:  # MTS law
                    taue = phase_data.get("taue", np.ones((nsyst, ngr_phase)))
                    for kkk in range(ngr_phase):
                        f.write(struct.pack(f"{nsyst}f", *taue[:, kkk]))

                # Twin data if present
                ntwmod = config.get(f"ntwmod_{iph}", 0)
                if ntwmod > 0:
                    twin_data = phase_data.get("twin_data", {})
                    pritw = twin_data.get("pritw", 0.0)
                    sectw = twin_data.get("sectw", 0.0)
                    ntwevents = twin_data.get(
                        "ntwevents", np.zeros(ngr_phase, dtype=int)
                    )
                    twfrph = twin_data.get("twfrph", np.zeros(ntwmod))
                    eftwfr = twin_data.get("eftwfr", np.zeros(ntwmod))
                    ktwsmx = twin_data.get("ktwsmx", np.zeros(ngr_phase, dtype=int))
                    twfrsy = twin_data.get("twfrsy", np.zeros((nsyst, ngr_phase)))

                    f.write(struct.pack("2f", pritw, sectw))
                    f.write(struct.pack(f"{ngr_phase}i", *ntwevents))
                    f.write(struct.pack(f"{ntwmod}f", *twfrph))
                    f.write(struct.pack(f"{ntwmod}f", *eftwfr))
                    f.write(struct.pack(f"{ngr_phase}i", *ktwsmx))

                    for kkk in range(ngr_phase):
                        f.write(struct.pack(f"{nsyst}f", *twfrsy[:, kkk]))

            print(f"Successfully wrote POSTMORT.OUT with full grain states")

    except (IOError, struct.error) as e:
        print(f"Error writing POSTMORT.OUT: {e}")
        raise
