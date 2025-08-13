#!/usr/bin/env python3
"""
Debug the ex03 process parsing
"""

input_path = "ex03_FCC/roll/vpsc8.in"

with open(input_path, "r") as f:
    lines = [
        line.strip() for line in f if line.strip() and not line.strip().startswith("*")
    ]

print("Non-comment lines:")
for i, line in enumerate(lines):
    print(f"{i:2d}: {line}")

print("\n" + "=" * 50)

# Find the line with number of processes
proc_idx = None
for idx, line in enumerate(lines):
    if line.isdigit() and idx > 20:  # Look after the main config section
        proc_idx = idx
        break

print(f"Found proc_idx at {proc_idx}: '{lines[proc_idx]}'")
if proc_idx:
    num_processes = int(lines[proc_idx])
    print(f"Number of processes: {num_processes}")

    # Parse each process block
    i = proc_idx + 1
    for proc_num in range(num_processes):
        print(f"\nProcess {proc_num+1}:")
        print(
            f"  Starting at line {i}: '{lines[i] if i < len(lines) else 'OUT OF BOUNDS'}'"
        )
        # Find ivgvar (first token of next non-empty, non-comment line)
        while i < len(lines) and (not lines[i] or lines[i].startswith("*")):
            i += 1
        if i >= len(lines):
            print(f"  Warning: Not enough lines for process {proc_num+1}, skipping.")
            continue
        try:
            ivgvar = int(lines[i].split()[0])
            print(f"  ivgvar = {ivgvar} from line: '{lines[i]}'")
        except Exception as e:
            print(f"  Error parsing ivgvar: {e} from line: '{lines[i]}'")
            i += 1
            continue
        i += 1

        # Find process file (next non-empty, non-comment, non-comment-only line)
        print(
            f"  Looking for process file starting at line {i}: '{lines[i] if i < len(lines) else 'OUT OF BOUNDS'}'"
        )
        while i < len(lines):
            tokens = lines[i].split()
            # Accept as process file if it looks like a filename (not just a number or comment)
            if tokens and not tokens[0].isdigit() and not lines[i].startswith("*"):
                proc_file = lines[i]
                print(f"  Found process file: '{proc_file}'")
                break
            print(f"  Skipping line {i}: '{lines[i]}'")
            i += 1
        else:
            print(f"  Warning: Could not find process file for process {proc_num+1}")
            continue
        i += 1

        # For ivgvar==2 or 3, find next non-empty, non-comment line for extra parameter
        if ivgvar == 2 and i < len(lines):
            print(
                f"  Looking for PCYS parameter at line {i}: '{lines[i] if i < len(lines) else 'OUT OF BOUNDS'}'"
            )
            while i < len(lines) and (not lines[i] or lines[i].startswith("*")):
                i += 1
            if i < len(lines):
                print(f"  Found PCYS parameter: '{lines[i]}'")
                i += 1
        if ivgvar == 3 and i < len(lines):
            print(
                f"  Looking for Lankford parameter at line {i}: '{lines[i] if i < len(lines) else 'OUT OF BOUNDS'}'"
            )
            while i < len(lines) and (not lines[i] or lines[i].startswith("*")):
                i += 1
            if i < len(lines):
                print(f"  Found Lankford parameter: '{lines[i]}'")
                i += 1
