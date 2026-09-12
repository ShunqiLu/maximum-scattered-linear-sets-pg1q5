from __future__ import annotations

import sympy as sp


def pfaffians_4_by_4(matrix: sp.Matrix) -> list[sp.Expr]:
    out: list[sp.Expr] = []
    for omitted in range(5):
        i, j, k, ell = [a for a in range(5) if a != omitted]
        out.append(
            sp.expand(
                matrix[i, j] * matrix[k, ell]
                - matrix[i, k] * matrix[j, ell]
                + matrix[i, ell] * matrix[j, k]
            )
        )
    return out


def wedge(u: sp.Matrix, v: sp.Matrix) -> sp.Matrix:
    out = sp.zeros(5)
    for i in range(5):
        for j in range(i + 1, 5):
            out[i, j] = sp.expand(u[i] * v[j] - u[j] * v[i])
            out[j, i] = -out[i, j]
    return out


def assert_zero(expr: sp.Expr, label: str) -> None:
    if sp.factor(sp.cancel(expr)) != 0:
        raise AssertionError(f"{label}: {sp.factor(sp.cancel(expr))}")


def assert_vector_zero(vector: sp.Matrix, label: str) -> None:
    for i, entry in enumerate(vector):
        assert_zero(entry, f"{label}[{i}]")


def assert_zero_mod_p(expr: sp.Expr, prime: int, label: str) -> None:
    numerator = sp.together(expr).as_numer_denom()[0]
    if sp.expand(numerator) == 0:
        return
    symbols = sorted(numerator.free_symbols, key=lambda item: item.name)
    if not symbols:
        if int(numerator) % prime != 0:
            raise AssertionError(f"{label} over GF({prime}): {numerator}")
        return
    polynomial = sp.Poly(sp.expand(numerator), *symbols, modulus=prime)
    if not polynomial.is_zero:
        raise AssertionError(
            f"{label} over GF({prime}): {polynomial.as_expr()}"
        )


def assert_nonzero_mod_p(value: int, prime: int, label: str) -> None:
    if value % prime == 0:
        raise AssertionError(f"{label} vanishes over GF({prime})")


def jacobian(quadrics: list[sp.Expr], variables: tuple[sp.Symbol, ...]) -> sp.Matrix:
    return sp.Matrix([[sp.diff(q, x) for x in variables] for q in quadrics])


def check_incidence_count() -> None:
    q = sp.symbols("Q", positive=True)
    p3 = (q**4 - 1) / (q - 1)
    p4 = (q**5 - 1) / (q - 1)
    grassmann = (q**5 - 1) * (q**4 - 1) / ((q**2 - 1) * (q - 1))
    i4 = p3 + q**2 * (q + 1) * (q**2 + 1)
    i2 = grassmann - q**6
    assert_zero(i2 - i4 - q**4, "incidence I2-I4")
    assert_zero(p4 * i4 - grassmann * p3, "incidence fibre identity")


def check_rank_one_trace_model() -> None:
    c = sp.symbols("c0:5", nonzero=True)
    t = sp.symbols("t0:5")

    alternating = sp.zeros(5)
    for i in range(5):
        nxt = (i + 1) % 5
        e_i = sp.eye(5).col(i)
        e_next_minus_trace = sp.eye(5).col(nxt) - c[i] * sp.ones(5, 1)
        alternating += t[i] * wedge(e_i, e_next_minus_trace)

    for i in range(5):
        assert_zero(alternating[i, i], f"trace matrix diagonal {i}")
    for i in range(5):
        for j in range(5):
            assert_zero(alternating[i, j] + alternating[j, i], f"trace alternating {i},{j}")

    c0, c1, c2, c3, c4 = c
    t0, t1, t2, t3, t4 = t
    stated_matrix = sp.Matrix([
        [0, (1-c0)*t0+c1*t1, -c0*t0+c2*t2, -c0*t0+c3*t3, -c0*t0+(c4-1)*t4],
        [(c0-1)*t0-c1*t1, 0, (1-c1)*t1+c2*t2, -c1*t1+c3*t3, -c1*t1+c4*t4],
        [c0*t0-c2*t2, (c1-1)*t1-c2*t2, 0, (1-c2)*t2+c3*t3, -c2*t2+c4*t4],
        [c0*t0-c3*t3, c1*t1-c3*t3, (c2-1)*t2-c3*t3, 0, (1-c3)*t3+c4*t4],
        [c0*t0+(1-c4)*t4, c1*t1-c4*t4, c2*t2-c4*t4, (c3-1)*t3-c4*t4, 0],
    ])
    for i in range(5):
        for j in range(5):
            assert_zero(alternating[i, j] - stated_matrix[i, j], f"trace displayed matrix {i},{j}")

    q = pfaffians_4_by_4(alternating)
    if len(q) != 5:
        raise AssertionError("trace model did not produce five principal Pfaffians")
    displayed = [
        -c1*t1*t2-c1*t1*t3+c2*t2*t3-c3*t1*t3+c4*t1*t4+c4*t2*t4+t1*t3,
        -c0*t0*t2-c0*t0*t3+c2*t2*t3+c2*t2*t4-c3*t3*t4+c4*t2*t4-t2*t4,
        -c0*t0*t3+c1*t1*t3+c1*t1*t4-c3*t0*t3-c3*t3*t4+c4*t0*t4+t0*t3,
        -c0*t0*t1+c1*t1*t4-c2*t0*t2-c2*t2*t4+c4*t0*t4+c4*t1*t4-t1*t4,
        -c0*t0*t1-c0*t0*t2+c1*t1*t2-c2*t0*t2+c3*t0*t3+c3*t1*t3+t0*t2,
    ]
    for i, (derived, stated) in enumerate(zip(q, displayed)):
        assert_zero(derived - stated, f"rank-one Pfaffian Q{i}")

    j = jacobian(q, t)
    at_e0 = {t[i]: int(i == 0) for i in range(5)}
    j0 = j.subs(at_e0)
    jacobian_minor = j0.extract([1, 2, 3], [1, 2, 4]).det()
    assert_zero(
        jacobian_minor - c0**2*c4,
        "trace Jacobian minor rows Q1,Q2,Q3 and columns t1,t2,t4",
    )

    exceptional = {c3: 1 - c0}
    v_exceptional = sp.Matrix([0, 1, -c0/c2, c0/c2, 0])
    assert_vector_zero((j0 * v_exceptional).subs(exceptional), "exceptional tangent")
    assert_zero(q[1].subs(dict(zip(t, v_exceptional))).subs(exceptional) + c0**2/c2,
                "exceptional Q1")

    d = c0 + c3 - 1
    b = c0 + c2 + c3 - 1
    v = sp.Matrix([0, c4*b/(c0*d), -c4/d, c4/d, 1])
    assert_vector_zero(j0 * v, "generic tangent")
    q_at_v = [sp.factor(poly.subs(dict(zip(t, v)))) for poly in q]
    assert_zero(q_at_v[4] + c4**2*(c1-c3)*b/(c0*d**2), "generic Q4")
    assert_zero(q_at_v[3].subs({c2: 1-c0-c3}) + c4, "generic boundary Q3")

    scalar = sp.symbols("c", nonzero=True)
    scalar_sub = {ci: scalar for ci in c}
    polynomial = 7*scalar**2 - 5*scalar + 1
    scalar_values = [sp.factor(value.subs(scalar_sub)) for value in q_at_v]
    expected = [
        0,
        -scalar*polynomial/(2*scalar-1)**2,
        scalar*polynomial/(2*scalar-1)**2,
        polynomial/(2*scalar-1),
        0,
    ]
    for i, (derived, stated) in enumerate(zip(scalar_values, expected)):
        assert_zero(derived - stated, f"scalar tangent Q{i}")

    d_scalar = 2*scalar - 1
    v0 = sp.Matrix([
        0,
        (3*scalar-1)/d_scalar,
        -scalar/d_scalar,
        scalar/d_scalar,
        1,
    ])
    v1 = sp.Matrix([v0[(i-1) % 5] for i in range(5)])
    e0 = sp.eye(5).col(0)
    e1 = sp.eye(5).col(1)
    line_matrix = sp.Matrix.hstack(e0, v0, e1, v1)
    consecutive_minor = line_matrix.extract([0, 1, 2, 3], range(4)).det()
    assert_zero(
        consecutive_minor - scalar/d_scalar,
        "consecutive Frobenius tangent-line minor",
    )

    delta_lp = scalar/(1-2*scalar)
    assert_zero(
        (1-2*scalar)**2*(delta_lp**2-delta_lp+1)-polynomial,
        "trace-root LP parameter",
    )
    assert_zero(
        (1-2*scalar)**2*(delta_lp*(1-delta_lp)-1)+polynomial,
        "trace-root LP determinant",
    )
    tau = sp.symbols("tau")
    trace_operator = sum(tau**i for i in range(5))
    f_trace = tau-scalar*trace_operator
    g_lp = tau**2+delta_lp*tau**3
    cleared = sp.together(
        -delta_lp*(f_trace+1)-g_lp*(f_trace+delta_lp)
    ).as_numer_denom()[0]
    remainder = sp.rem(cleared, tau**5-1, tau)
    assert_zero(remainder+tau**3*polynomial, "trace-root LP operator identity")


def check_vertical_model() -> None:
    n = sp.symbols("n0:5", nonzero=True)
    z = sp.symbols("z0:5")
    standard = [sp.eye(5).col(i) for i in range(5)]
    r = standard[:4] + [-sum(standard[:4], sp.zeros(5, 1))]
    s = standard[4]

    alternating = sp.zeros(5)
    for i in range(5):
        alternating += z[i] * wedge(r[i], r[(i-1) % 5] + n[i]*s)
    for i in range(5):
        assert_zero(alternating[i, i], f"vertical matrix diagonal {i}")
    for i in range(5):
        for j in range(5):
            assert_zero(alternating[i, j] + alternating[j, i], f"vertical alternating {i},{j}")

    n0, n1, n2, n3, n4 = n
    z0, z1, z2, z3, z4 = z
    stated_matrix = sp.Matrix([
        [0, -z0-z1, -z0, -z0-z4, n0*z0-n4*z4],
        [z0+z1, 0, -z2, -z4, n1*z1-n4*z4],
        [z0, z2, 0, -z3-z4, n2*z2-n4*z4],
        [z0+z4, z4, z3+z4, 0, n3*z3-n4*z4],
        [-n0*z0+n4*z4, -n1*z1+n4*z4, -n2*z2+n4*z4, -n3*z3+n4*z4, 0],
    ])
    for i in range(5):
        for j in range(5):
            assert_zero(alternating[i, j] - stated_matrix[i, j], f"vertical displayed matrix {i},{j}")

    q = pfaffians_4_by_4(alternating)
    if len(q) != 5:
        raise AssertionError("vertical model did not produce five principal Pfaffians")

    displayed = [
        -n1*z1*z3-n1*z1*z4+n2*z2*z4-n3*z2*z3+n4*z2*z4+n4*z3*z4,
        -n0*z0*z3-n0*z0*z4+n2*z0*z2+n2*z2*z4-n3*z0*z3+n4*z3*z4,
        -n0*z0*z4+n1*z0*z1+n1*z1*z4-n3*z0*z3-n3*z1*z3+n4*z1*z4,
        -n0*z0*z2+n1*z0*z1-n2*z0*z2-n2*z1*z2+n4*z1*z4+n4*z2*z4,
        z0*z2+z0*z3+z1*z3+z1*z4+z2*z4,
    ]
    for i, (derived, stated) in enumerate(zip(q, displayed)):
        assert_zero(derived - stated, f"vertical Pfaffian Q{i}")

    j = jacobian(q, z)
    at_e0 = {z[i]: int(i == 0) for i in range(5)}
    j0 = j.subs(at_e0)
    jacobian_minor = j0.extract([1, 2, 3], [1, 2, 4]).det()
    assert_zero(
        jacobian_minor - n0**2*n1,
        "vertical Jacobian minor rows Q1,Q2,Q3 and columns z1,z2,z4",
    )

    delta = n0 + n2 + n3
    v = sp.Matrix([0, (n0+n2)/n1, 1, -1, delta/n0])
    trace_relation = {n4: -(n0+n1+n2+n3)}
    assert_vector_zero((j0 * v).subs(trace_relation), "vertical tangent")
    assert_zero(
        (q[1].subs(dict(zip(z, v))) - (n2-n4)*delta/n0).subs(trace_relation),
        "vertical Q1",
    )
    scalar_v = {z0: 0, z1: 2, z2: 1, z3: -1, z4: 3}
    assert_zero(q[4].subs(scalar_v) - 7, "vertical characteristic-five Q4")


def check_c4_scalar_identity() -> None:
    k, x = sp.symbols("k x")
    p = 7*x**2 - 5*x + 1

    tau = sp.symbols("tau")
    trace_operator = sum(tau**i for i in range(5))
    f_k = k*(tau + tau**3) + tau**2 + tau**4
    assert_zero((tau + k)*(1 + tau**2)*tau - f_k, "C4 Ore factorization")

    odd_inverse = trace_operator/2 - 1 - tau
    inverse_residual = sp.rem(
        sp.together((tau + tau**3)*odd_inverse - 1),
        tau**5 - 1,
        tau,
    )
    assert_zero(inverse_residual, "C4 odd inverse operator")

    assert_zero(
        (2*k+3)**2 * p.subs(x, (k+1)/(2*k+3)) - (k**2+k+1),
        "C4 scalar identity",
    )
    assert_zero(32*((-sp.Rational(3, 2))**5 - 1) + 275, "C4 char-11 numerator")
    assert_zero(25*p.subs(x, sp.Rational(1, 5)) - 7, "C4 P(1/5)")


def check_special_characteristics() -> None:
    x, k = sp.symbols("x k")
    polynomial = 7*x**2 - 5*x + 1

    assert_zero_mod_p((2*x - 1) - 1, 2, "trace denominator is a unit")
    assert_zero_mod_p(
        polynomial - (x**2 + x + 1),
        2,
        "characteristic-two polynomial",
    )
    assert_zero_mod_p(
        x*(1+x) - 1 - (x**2+x+1),
        2,
        "characteristic-two LP norm relation",
    )
    assert_zero_mod_p(
        (x**2+x+1).subs(x, k+1) - (k**2+k+1),
        2,
        "characteristic-two C4 scalar identity",
    )
    assert_nonzero_mod_p(7-5+1, 2, "characteristic-two P(1)")

    assert_zero_mod_p(
        polynomial - (x-1)**2,
        3,
        "characteristic-three polynomial",
    )

    assert_zero_mod_p(
        polynomial - (2*x**2+1),
        5,
        "characteristic-five polynomial",
    )
    assert_zero_mod_p(k**5-1-(k-1)**5, 5, "characteristic-five fifth powers")
    assert_nonzero_mod_p(7, 5, "characteristic-five vertical residual")

    assert_zero_mod_p(
        polynomial - 2*(x-3),
        7,
        "characteristic-seven polynomial",
    )

    k11 = (-3*pow(2, -1, 11)) % 11
    c11 = pow(5, -1, 11)
    if pow(k11, 5, 11) != 1:
        raise AssertionError("characteristic-eleven fifth-power specialization")
    assert_nonzero_mod_p(7*c11**2-5*c11+1, 11, "characteristic-eleven P(1/5)")


def run_symbolic_checks() -> list[tuple[str, str]]:
    checks: list[tuple[str, str]] = []
    check_incidence_count()
    checks.append(("INC", "rank-two/rank-four incidence identities"))
    check_rank_one_trace_model()
    checks.extend([
        ("TR-MAT/PF", "trace matrix and principal Pfaffians"),
        ("TR-JAC/TAN", "trace Jacobian and tangent evaluations"),
        ("TR-LINE/LP", "trace line minor and LP identity"),
    ])
    check_vertical_model()
    checks.extend([
        ("VE-MAT/PF", "vertical matrix and principal Pfaffians"),
        ("VE-JAC/TAN", "vertical Jacobian and tangent evaluations"),
    ])
    check_c4_scalar_identity()
    checks.append(("C4-OP/SCALAR", "C4 operator and scalar identities, including the two exceptional scalar calculations from Appendix B"))
    check_special_characteristics()
    checks.append(("SPECIAL-CHAR", "terminal specializations in characteristics 2, 3, 5, 7, and 11"))
    return checks


def main() -> None:
    for label, description in run_symbolic_checks():
        print(f"PASS [{label}] {description}")


if __name__ == "__main__":
    main()
