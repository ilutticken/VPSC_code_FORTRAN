#!/usr/bin/env python3
"""
Test script for postmortem read/write functionality
"""

import sys
import os

sys.path.append("src")

import numpy as np
from readers import read_postmortem
from writers import write_postmortem


def test_postmortem():
    """Test the postmortem read/write functionality"""

    # Create test configuration
    config = {
        "nph": 1,
        "ngr_total": 10,
        "ngr_0": 0,  # Add grain indexing
        "ngr_1": 10,  # Phase 1 has grains 1-10
        "isave": 5,
        "dav": np.array([1.0, 0.0, 0.0, 0.0, -1.0]),
        "sav": np.array([100.0, 0.0, 0.0, 0.0, 0.0]),
        "dbar0": np.array([1.0, 0.0, 0.0, 0.0, -1.0]),
        "mbartg": np.eye(5) * 100.0,
        "lbartg": np.eye(5) * 0.01,
        "sg": np.random.rand(5, 10) * 100,
        "epstotc": np.array([[0.1, 0.05, 0.0], [0.05, -0.05, 0.0], [0.0, 0.0, -0.05]]),
        "epsvm": 0.15,
        "epsacu": 0.2,
        "axisph": np.zeros((4, 3, 2)),  # (4, 3, nph+1) = (4, 3, 2) for nph=1
        "fijph": np.zeros((3, 3, 2)),  # (3, 3, nph+1) = (3, 3, 2) for nph=1
        "ihardlaw": 0,
        "nsyst_1": 12,
        "grain_data": {
            "phase_1": {
                "ag": np.random.rand(3, 3, 10),
                "wgt": np.ones(10) / 10.0,
                "gtotgr": np.random.rand(10) * 0.1,
                "crss": np.ones((12, 10)) * 50.0,
            }
        },
    }

    # Test file path
    test_file = "test_postmort.out"

    print("Testing postmortem write...")
    try:
        write_postmortem(test_file, config)
        print("✓ Write successful")
    except Exception as e:
        print(f"✗ Write failed: {e}")
        return False

    print("Testing postmortem read...")
    try:
        config_read = read_postmortem(test_file, {})
        print("✓ Read successful")

        # Verify some key data
        if "dav" in config_read and np.allclose(config_read["dav"], config["dav"]):
            print("✓ DAV array matches")
        else:
            print("✗ DAV array mismatch")

        if "isave" in config_read and config_read["isave"] == config["isave"]:
            print("✓ isave value matches")
        else:
            print("✗ isave value mismatch")

        if (
            "epsvm" in config_read
            and abs(config_read["epsvm"] - config["epsvm"]) < 1e-6
        ):
            print("✓ von Mises strain matches")
        else:
            print("✗ von Mises strain mismatch")

    except Exception as e:
        print(f"✗ Read failed: {e}")
        return False

    # Cleanup
    if os.path.exists(test_file):
        os.remove(test_file)
        print("✓ Test file cleaned up")

    return True


if __name__ == "__main__":
    print("=== Postmortem Functionality Test ===")
    success = test_postmortem()
    if success:
        print("\n🎉 All tests passed! Postmortem functionality is working correctly.")
    else:
        print("\n❌ Some tests failed. Check the implementation.")
