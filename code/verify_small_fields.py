from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

sys.dont_write_bytecode = True

import galois
import numpy as np


Q_VALUES = (2, 3, 4, 5, 7, 8, 9)
EXHAUSTIVE_Q = frozenset((2, 3, 4, 5))
SAMPLE_SIZE = 64
LLZ_Q_VALUES = (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25)


@dataclass(frozen=True)
class FieldCheck:
    q: int
    trace_tested: int
    trace_total: int
    trace_scattered: int
    trace_expected_scattered: int
    c3_tested: int
    c3_total: int
    c4_tested: int
    c4_total: int
    exhaustive: bool


@dataclass(frozen=True)
class LLZRangeCheck:
    q: int
    c3_excluded: int
    c3_total: int
    c3_shift_cutoff: int
    c4_excluded: int
    c4_total: int
    c4_shift_cutoff: int


def relative_trace(field: type, q: int, values):
    out = field.Zeros(np.shape(values))
    conjugate = values
    for _ in range(5):
        out += conjugate
        conjugate = conjugate**q
    return out


def projective_fibre_counts(first, second) -> np.ndarray:
    nonzero = (first != 0) | (second != 0)
    if int(np.count_nonzero(nonzero)) != first.size - 1:
        raise AssertionError("the displayed parameter map is not injective")
    x = first[nonzero]
    y = second[nonzero]
    finite = x != 0
    codes = np.empty(x.size, dtype=np.int64)
    codes[~np.asarray(finite)] = 0
    codes[np.asarray(finite)] = 1 + np.asarray(y[finite] / x[finite], dtype=np.int64)
    return np.unique(codes, return_counts=True)[1]


def graph_is_scattered(elements, values, q: int) -> bool:
    counts = projective_fibre_counts(elements, values)
    return bool(np.all(counts == q - 1))


def subspace_is_scattered(first, second, q: int) -> bool:
    counts = projective_fibre_counts(first, second)
    return bool(np.all(counts == q - 1))


def deterministic_sample(values, size: int):
    if values.size <= size:
        return values
    indices = np.linspace(0, values.size - 1, num=size, dtype=np.int64)
    return values[np.unique(indices)]


def union_parameters(field: type, *arrays):
    encoded = np.concatenate([np.asarray(array, dtype=np.int64) for array in arrays])
    return field(np.unique(encoded))


def trace_prediction(c, q: int) -> bool:
    field = type(c)
    characteristic = field.characteristic
    polynomial = (
        field(7 % characteristic)*c**2
        - field(5 % characteristic)*c
        + field(1)
    )
    return bool(c == 0 or (c**q == c and c != 0 and polynomial == 0))


def check_one_field(q: int) -> FieldCheck:
    field = galois.GF(q**5)
    elements = field.elements
    trace_values = relative_trace(field, q, elements)
    scalar_mask = elements**q == elements
    scalars = elements[scalar_mask]

    if q in EXHAUSTIVE_Q:
        trace_parameters = elements
    else:
        nonscalars = elements[~scalar_mask]
        trace_parameters = union_parameters(
            field,
            scalars,
            deterministic_sample(nonscalars, SAMPLE_SIZE),
        )

    x_q = elements**q
    observed_scattered = 0
    expected_scattered = 0
    for c in trace_parameters:
        expected = trace_prediction(c, q)
        actual = graph_is_scattered(elements, x_q - c*trace_values, q)
        if actual != expected:
            raise AssertionError(
                f"q={q}: trace parameter {int(c)} has actual={actual}, expected={expected}"
            )
        observed_scattered += int(actual)
        expected_scattered += int(expected)

    trace_zero_nonzero = elements[(trace_values == 0) & (elements != 0)]
    c3_parameters = (
        trace_zero_nonzero
        if q in EXHAUSTIVE_Q
        else deterministic_sample(trace_zero_nonzero, SAMPLE_SIZE)
    )
    rho = next(value for value in elements if relative_trace(field, q, value) != 0)
    t_values = relative_trace(field, q, rho*elements)
    y_values = x_q - elements
    second_c3 = x_q - elements ** (q**4)
    for eta in c3_parameters:
        first_c3 = eta*y_values + t_values
        if subspace_is_scattered(first_c3, second_c3, q):
            raise AssertionError(f"q={q}: C3 parameter eta={int(eta)} is scattered")

    norm_exponent = (q**5 - 1) // (q - 1)
    c4_all = elements[(elements != 0) & (elements**norm_exponent == 1)]
    c4_parameters = c4_all if q in EXHAUSTIVE_Q else deterministic_sample(c4_all, SAMPLE_SIZE)
    x_q2 = elements ** (q**2)
    x_q3 = elements ** (q**3)
    x_q4 = elements ** (q**4)
    for k in c4_parameters:
        values = k*(x_q + x_q3) + x_q2 + x_q4
        if graph_is_scattered(elements, values, q):
            raise AssertionError(f"q={q}: C4 parameter k={int(k)} is scattered")

    return FieldCheck(
        q=q,
        trace_tested=int(trace_parameters.size),
        trace_total=int(elements.size),
        trace_scattered=observed_scattered,
        trace_expected_scattered=expected_scattered,
        c3_tested=int(c3_parameters.size),
        c3_total=q**4 - 1,
        c4_tested=int(c4_parameters.size),
        c4_total=int(c4_all.size),
        exhaustive=q in EXHAUSTIVE_Q,
    )


def run_small_field_checks() -> list[FieldCheck]:
    return [check_one_field(q) for q in Q_VALUES]


def nonzero_trace_element(field: type, q: int):
    rho = field(1)
    if bool(relative_trace(field, q, rho) != 0):
        return rho
    alpha = field.primitive_element
    rho = alpha
    while bool(relative_trace(field, q, rho) == 0):
        rho = rho*alpha
    return rho


def c3_collision_certificate(field: type, q: int) -> tuple[int, int, int]:
    alpha = field.primitive_element
    quotient_order = (q**5-1)//(q-1)
    representatives = alpha**np.arange(quotient_order, dtype=np.int64)
    rho = nonzero_trace_element(field, q)
    first_linear = representatives**q-representatives
    second = representatives**q-representatives**(q**4)
    trace_linear = relative_trace(field, q, rho*representatives)
    excluded: set[int] = set()
    total = q**4-1
    for shift in range(1, quotient_order):
        first_shift = np.roll(first_linear, -shift)
        second_shift = np.roll(second, -shift)
        trace_shift = np.roll(trace_linear, -shift)
        denominator = second*first_shift-second_shift*first_linear
        numerator = second_shift*trace_linear-second*trace_shift
        zero_denominator = denominator == 0
        if bool(np.any(zero_denominator & (numerator == 0))):
            return total, total, shift
        valid = ~zero_denominator
        eta = numerator[valid]/denominator[valid]
        admissible = (eta != 0) & (relative_trace(field, q, eta) == 0)
        excluded.update(np.asarray(eta[admissible], dtype=np.int64).tolist())
        if len(excluded) == total:
            return len(excluded), total, shift
        if shift % 100 == 0:
            print(
                f"PROGRESS [LLZ-C3] q={q} shift={shift} "
                f"coverage={len(excluded)}/{total}",
                file=sys.stderr,
                flush=True,
            )
    raise AssertionError(f"q={q}: C3 collision coverage {len(excluded)}/{total}")


def c4_collision_certificate(field: type, q: int) -> tuple[int, int, int]:
    alpha = field.primitive_element
    total = (q**5-1)//(q-1)
    y = alpha**((q-1)*np.arange(total, dtype=np.int64))
    coefficient = y+y**(q**2+q+1)
    constant = y**(q+1)+y**(q**3+q**2+q+1)
    excluded: set[int] = set()
    for shift in range(1, total):
        coefficient_difference = np.roll(coefficient, -shift)-coefficient
        constant_difference = np.roll(constant, -shift)-constant
        zero_denominator = coefficient_difference == 0
        if bool(np.any(zero_denominator & (constant_difference == 0))):
            return total, total, shift
        valid = ~zero_denominator
        parameter = -constant_difference[valid]/coefficient_difference[valid]
        parameter = parameter[parameter != 0]
        logarithms = parameter.log()
        norm_one = logarithms % (q-1) == 0
        indices = logarithms[norm_one]//(q-1)
        excluded.update(np.asarray(indices, dtype=np.int64).tolist())
        if len(excluded) == total:
            return len(excluded), total, shift
        if shift % 100 == 0:
            print(
                f"PROGRESS [LLZ-C4] q={q} shift={shift} "
                f"coverage={len(excluded)}/{total}",
                file=sys.stderr,
                flush=True,
            )
    raise AssertionError(f"q={q}: C4 collision coverage {len(excluded)}/{total}")


def check_llz_field(q: int) -> LLZRangeCheck:
    print(f"PROGRESS [LLZ] q={q} constructing GF({q}^5)", file=sys.stderr, flush=True)
    field = galois.GF(q**5, compile="jit-lookup")
    c3_excluded, c3_total, c3_shift = c3_collision_certificate(field, q)
    print(
        f"PROGRESS [LLZ-C3] q={q} complete at shift {c3_shift}",
        file=sys.stderr,
        flush=True,
    )
    c4_excluded, c4_total, c4_shift = c4_collision_certificate(field, q)
    print(
        f"PROGRESS [LLZ-C4] q={q} complete at shift {c4_shift}",
        file=sys.stderr,
        flush=True,
    )
    return LLZRangeCheck(
        q=q,
        c3_excluded=c3_excluded,
        c3_total=c3_total,
        c3_shift_cutoff=c3_shift,
        c4_excluded=c4_excluded,
        c4_total=c4_total,
        c4_shift_cutoff=c4_shift,
    )


def run_llz_range_checks() -> list[LLZRangeCheck]:
    script = Path(__file__).resolve()
    results: list[LLZRangeCheck] = []
    for q in LLZ_Q_VALUES:
        print(f"RUN [LLZ-range] q={q}", flush=True)
        completed = subprocess.run(
            [sys.executable, "-B", str(script), "--llz-worker", str(q)],
            stdout=subprocess.PIPE,
            text=True,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(f"q={q}: LLZ worker failed with code {completed.returncode}")
        results.append(LLZRangeCheck(**json.loads(completed.stdout)))
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--llz-worker", type=int)
    args = parser.parse_args()
    if args.llz_worker is not None:
        if args.llz_worker not in LLZ_Q_VALUES:
            raise ValueError(f"unsupported q={args.llz_worker}")
        print(json.dumps(asdict(check_llz_field(args.llz_worker))))
        return
    for result in run_small_field_checks():
        scope = "exhaustive" if result.exhaustive else "deterministic sample"
        print(
            "PASS [small-fields] "
            f"q={result.q} ({scope}): trace {result.trace_tested}/{result.trace_total}, "
            f"C3 {result.c3_tested}/{result.c3_total}, "
            f"C4 {result.c4_tested}/{result.c4_total}"
        )
    for result in run_llz_range_checks():
        print(
            "PASS [LLZ-range] "
            f"q={result.q}: C3 {result.c3_excluded}/{result.c3_total} "
            f"by shift {result.c3_shift_cutoff}, "
            f"C4 {result.c4_excluded}/{result.c4_total} "
            f"by shift {result.c4_shift_cutoff}"
        )


if __name__ == "__main__":
    main()
