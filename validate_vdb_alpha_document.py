#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Validate the generalized VdB + generic-lapse formulas against GR_python.

The metric under test is

    ds^2 = -alpha(r)^2 dt^2
           + B(r)^2[(dr - beta(r) dt)^2 + r^2 dOmega^2].

This script keeps GR_python's orthonormal-frame convention for the radial flux:

    8*pi*j_hat{r} = -(2*beta/(alpha*B))*V_tilde.

Usage:
    python validate_vdb_alpha_document.py
"""

import sympy as sp
from sympy import Function, Matrix, diff, pi, sin, symbols

import gr_main as gm
from gr_metric_library import select_metric
from gr_warp import build_vdb_alpha_metric_configuration, document_vdb_alpha_formulas


def simplify_residual(expr):
    if isinstance(expr, sp.MatrixBase):
        return Matrix(expr.rows, expr.cols, lambda i, j: simplify_residual(expr[i, j]))
    return sp.simplify(sp.expand(sp.cancel(sp.together(expr))))


def is_zero(expr):
    if isinstance(expr, sp.MatrixBase):
        return all(is_zero(expr[i, j]) for i in range(expr.rows) for j in range(expr.cols))
    return simplify_residual(expr) == 0


def check(label, computed, expected):
    try:
        residual = simplify_residual(computed - expected)
    except TypeError:
        residual = 0 if computed == expected else f"{computed!r} != {expected!r}"
    passed = is_zero(residual)
    print(f"[{'MATCH' if passed else 'MISMATCH'}] {label}")
    print(f"  computed : {computed}")
    print(f"  expected : {expected}")
    print(f"  residual : {residual}")
    return {
        "label": label,
        "computed": computed,
        "expected": expected,
        "residual": residual,
        "passed": passed,
    }


def check_structural(label, computed, expected):
    computed_render = sp.sstr(computed)
    expected_render = sp.sstr(expected)
    passed = computed_render == expected_render
    residual = 0 if passed else "render mismatch"
    print(f"[{'MATCH' if passed else 'MISMATCH'}] {label}")
    print(f"  computed : {computed_render}")
    print(f"  expected : {expected_render}")
    print(f"  residual : {residual}")
    return {
        "label": label,
        "computed": computed_render,
        "expected": expected_render,
        "residual": residual,
        "passed": passed,
    }


def run_context():
    t, r, theta, phi = symbols("t r theta phi", real=True)
    coords = [t, r, theta, phi]
    alpha = Function("alpha", positive=True)(r)
    B = Function("B", positive=True)(r)
    beta = Function("beta")(r)
    v = symbols("v", real=True)

    cfg = build_vdb_alpha_metric_configuration(coords, alpha, beta, B)
    formulas = document_vdb_alpha_formulas(r, alpha, B, beta, v)
    results = gm.run_computations(
        g_metric=cfg["g_metric"],
        coords=coords,
        dim=4,
        g_inv_metric=None,
        e_tetrad=None,
        compute_weyl_flag=False,
        compute_kretschmann_flag=False,
        compute_geodesics_flag=False,
        compute_killing_flag=False,
        compute_tetrad_flag=True,
        fast_mode=True,
        compute_horizons_flag=False,
        compute_penrose_flag=False,
    )
    return r, theta, alpha, B, beta, v, cfg, formulas, results


def collect_checks():
    r, theta, alpha, B, beta, v, cfg, formulas, results = run_context()
    builtin = select_metric("warp_doc_variant_b_alpha", [symbols("t", real=True), r, theta, symbols("phi", real=True)])

    G_ortho = results["G_ortho"]
    rho = sp.cancel(G_ortho[0, 0] / (8 * pi))
    j_hat = sp.cancel(G_ortho[0, 1] / (8 * pi))
    p_r = sp.cancel(G_ortho[1, 1] / (8 * pi))
    rho_plus_pr = sp.cancel(rho + p_r)

    alpha_magic = formulas["alpha_magic"]
    alpha_magic_subs = {
        alpha: alpha_magic,
        diff(alpha, r): diff(alpha_magic, r),
        diff(alpha, r, 2): diff(alpha_magic, r, 2),
    }

    gamma_sq = sp.cancel(1 / (1 - v**2))
    boosted_sum = sp.cancel(gamma_sq * ((1 + v**2) * rho_plus_pr - 4 * v * j_hat))

    checks = [
        check_structural("built-in metric registration", builtin["g_metric"], cfg["g_metric"]),
        check_structural("built-in tetrad registration", builtin["e_tetrad"], cfg["e_tetrad"]),
        check("rho", rho, formulas["rho"]),
        check("j_hat{r} with GR_python sign", j_hat, formulas["j_r"]),
        check("rho + p_r", rho_plus_pr, formulas["rho_plus_pr"]),
        check("boosted rho_plus_pr factorization", boosted_sum, formulas["boosted_sum_factorized"]),
        check("alpha_magic = (rB)'/B", alpha_magic, sp.cancel(diff(r * B, r) / B)),
        check("V_tilde(alpha_magic)", formulas["V_tilde"].subs(alpha_magic_subs), 0),
        check("j_hat{r}(alpha_magic)", j_hat.subs(alpha_magic_subs), 0),
        check("(rho + p_r)(alpha_magic)", rho_plus_pr.subs(alpha_magic_subs), 0),
        check("first-integral differential form", formulas["first_integral_lhs"], formulas["first_integral_rhs"]),
    ]
    return checks


def main():
    checks = collect_checks()
    passed = sum(item["passed"] for item in checks)
    total = len(checks)
    print(f"\nChecks passed: {passed}/{total}")
    raise SystemExit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
