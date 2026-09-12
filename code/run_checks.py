"""Run the mathematical checks and compare their results with the reference data."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import platform
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True

import galois
import llvmlite
import numba
import numpy as np
import sympy as sp

from verify_small_fields import (
    LLZ_Q_VALUES,
    Q_VALUES,
    LLZRangeCheck,
    run_small_field_checks,
)
from verify_symbolic_identities import run_symbolic_checks


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_VERSIONS = {
    "python": "3.12.13",
    "sympy": "1.14.0",
    "galois": "0.4.11",
    "numpy": "2.3.5",
    "numba": "0.67.0",
    "llvmlite": "0.49.0",
}


def check_environment() -> None:
    actual = {
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "galois": galois.__version__,
        "numpy": np.__version__,
        "numba": numba.__version__,
        "llvmlite": llvmlite.__version__,
    }
    if actual != EXPECTED_VERSIONS:
        differences = ", ".join(
            f"{name}: expected {EXPECTED_VERSIONS[name]}, got {actual[name]}"
            for name in EXPECTED_VERSIONS
            if actual[name] != EXPECTED_VERSIONS[name]
        )
        raise RuntimeError(f"verification environment mismatch: {differences}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run all symbolic and small-field checks, with LLZ collisions limited to q <= 9.",
    )
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError("Run without -O or PYTHONOPTIMIZE so all checks remain enabled.")
    check_environment()
    print("PASS [environment] locked interpreter and dependency versions", flush=True)
    reference = json.loads((ROOT / "data/reference_results.json").read_text(encoding="utf-8"))
    if [row["q"] for row in reference["small_fields"]] != list(Q_VALUES):
        raise AssertionError("reference small-field coverage differs")
    if [row["q"] for row in reference["llz_fields"]] != list(LLZ_Q_VALUES):
        raise AssertionError("reference LLZ coverage differs")

    symbolic = run_symbolic_checks()
    if [label for label, _ in symbolic] != reference["symbolic_check_ids"]:
        raise AssertionError("symbolic check groups differ from reference")
    for label, description in symbolic:
        print(f"PASS [{label}] {description}", flush=True)

    fields = run_small_field_checks()
    if [asdict(result) for result in fields] != reference["small_fields"]:
        raise AssertionError("small-field results differ from reference")
    for result in fields:
        print(f"PASS [small-fields] q={result.q}: exact reference match", flush=True)

    selected = [q for q in LLZ_Q_VALUES if not args.quick or q <= 9]
    script = Path(__file__).with_name("verify_small_fields.py")
    expected_llz = {row["q"]: row for row in reference["llz_fields"]}
    for q in selected:
        print(f"RUN [LLZ-range] q={q}", flush=True)
        completed = subprocess.run(
            [sys.executable, "-B", str(script), "--llz-worker", str(q)],
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )
        result = LLZRangeCheck(**json.loads(completed.stdout))
        if asdict(result) != expected_llz[q]:
            raise AssertionError(f"q={q}: collision results differ from reference")
        print(f"PASS [LLZ-range] q={q}: exact reference match", flush=True)

    print(
        f"PASS: {len(symbolic)} symbolic groups, {len(fields)} small-field rows, "
        f"{len(selected)}/{len(LLZ_Q_VALUES)} LLZ collision rows.",
        flush=True,
    )
    if args.quick:
        print("Quick mode completed; run without --quick to check every prime power q <= 25.")


if __name__ == "__main__":
    main()
