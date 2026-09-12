# Reproducibility data

`reference_results.json` contains the expected results for
[`../code/run_checks.py`](../code/run_checks.py):

- `symbolic_check_ids`: the eight symbolic groups that must pass;
- `small_fields`: seven finite-field rows with parameter counts, observed
  and predicted scattered trace counts, and the exhaustive or sampled scope;
- `llz_fields`: fourteen collision-certificate rows, one for each prime
  power `q <= 25`, recording excluded C3 and C4 parameters and the cyclic
  shift at which each enumeration achieves full coverage.

The values are the reference numerical results accompanying Version 1 of
the preprint. All numerical comparisons are exact. The entry point reads
this file and recomputes the selected checks without modifying it.
