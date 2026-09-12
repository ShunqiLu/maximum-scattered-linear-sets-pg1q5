# Reproducibility code

The two mathematical modules verify symbolic identities and finite-field
collision counts. They use the public reference data in
[`../data/reference_results.json`](../data/reference_results.json).

## Environment and command

`pyproject.toml` requires CPython 3.12.13 and pins the direct dependencies.
`uv.lock` records the complete resolved dependency graph. With
[uv](https://docs.astral.sh/uv/) installed, run from the repository root:

```text
uv run --project code --locked python -B code/run_checks.py
```

The same entry point also works from this directory:

```text
uv run --locked python -B run_checks.py
```

The full check can take substantial time because it constructs collision
certificates for every admissible C3 and C4 parameter for every prime power
`q <= 25`. An optional shorter run checks the same eight symbolic groups
and all seven small-field rows, with collision enumeration limited to `q <= 9`:

```text
uv run --project code --locked python -B code/run_checks.py --quick
```

The programs print progress and compare the integer results exactly with
the reference data. A failed identity, result mismatch or environment
mismatch exits with an error. No input files beyond the public reference
data are required.

## Calculations

- `verify_symbolic_identities.py`: regenerates the two alternating matrices
  and their principal Pfaffians, then checks the Jacobian minors, tangent
  substitutions, line minor, operator identities and specializations in
  characteristics 2, 3, 5, 7 and 11.
- `verify_small_fields.py`: checks the trace, C3 and C4 families over small
  finite fields and constructs direct collision certificates throughout
  the LLZ range.
- `run_checks.py`: checks the environment, runs the mathematical modules
  and compares the results with the reference data.

For `q = 2, 3, 4, 5`, the small-field checks exhaust all parameters. For
`q = 7, 8, 9`, they test every scalar trace parameter and a fixed evenly
spaced sample of the other parameters. The separate LLZ calculation
exhausts every admissible C3 and C4 parameter for each prime power
`q <= 25`. These finite computations supplement the proof; the general
theorems do not rely on enumeration.
