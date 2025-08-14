#!/usr/bin/env python3
"""
VPSC8 Path Fixer - Convert Windows-style paths to cross-platform paths in input files.

This script fixes path separators in VPSC8 input files to ensure cross-platform compatibility.
Run this script if you encounter path-related issues on Linux/macOS systems.

Usage:
    python fix_paths.py [directory]

If no directory is specified, it will process all example directories.
"""

import os
import glob
import sys


def fix_paths_in_file(filepath):
    """Fix Windows-style paths in a single file."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()

        # Replace Windows-style paths with cross-platform paths
        # Pattern: ex##_name\file.ext -> file.ext (just use local filename)
        import re

        # Find patterns like ex02_FCC\file.ext and replace with just file.ext
        pattern = r"ex\d+_\w+\\([^\\]+)"
        content = re.sub(pattern, r"\1", content)

        # Write back if changes were made
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(content)

        return True
    except Exception as e:
        print(f"Warning: Could not process {filepath}: {e}")
        return False


def fix_paths_in_directory(directory):
    """Fix paths in all .in files in a directory."""
    input_files = glob.glob(os.path.join(directory, "*.in"))
    input_files.extend(glob.glob(os.path.join(directory, "*", "*.in")))

    fixed_count = 0
    for filepath in input_files:
        if fix_paths_in_file(filepath):
            fixed_count += 1
            print(f"Fixed: {filepath}")

    return fixed_count


def main():
    if len(sys.argv) > 1:
        directory = sys.argv[1]
        if not os.path.exists(directory):
            print(f"Error: Directory {directory} does not exist")
            sys.exit(1)
        directories = [directory]
    else:
        # Process all example directories
        examples_dir = "examples"
        if not os.path.exists(examples_dir):
            print(
                "Error: examples directory not found. Run this script from the VPSC8 root directory."
            )
            sys.exit(1)

        directories = [
            os.path.join(examples_dir, d)
            for d in os.listdir(examples_dir)
            if os.path.isdir(os.path.join(examples_dir, d)) and d.startswith("ex")
        ]

    total_fixed = 0
    for directory in directories:
        print(f"Processing directory: {directory}")
        fixed = fix_paths_in_directory(directory)
        total_fixed += fixed

    print(f"\nCompleted! Fixed {total_fixed} files.")
    print("All input files now use local filenames for cross-platform compatibility.")


if __name__ == "__main__":
    main()
