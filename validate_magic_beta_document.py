#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Validate the later "Magic Beta" identities against GR_python.

This script assumes the generic-lapse VdB metric and then imposes the magic
lapse alpha = (rB)'/B. It verifies the algebraic identities behind the
positivity-barrier and causal-structure discussion.

Usage:
    python validate_magic_beta_document.py
"""

import sympy as sp
from sympy import diff, pi, sin

from validate_vdb_alpha_document import check, run_context


def collect_checks():
    r, theta, alpha, B, beta, v, cfg, formulas, results = run_context()

    g = cfg["g_metric"]
    g_inv = sp.simplify(g.inv())
    gamma = sp.diag(B**2, B**2 * r**2, B**2 * r**2 * sin(theta) ** 2)

    G_ortho = results["G_ortho"]
    rho = sp.cancel(G_ortho[0, 0] / (8 * pi))
    j_hat = sp.cancel(G_ortho[0, 1] / (8 * pi))
    p_r = sp.cancel(G_ortho[1, 1] / (8 * pi))

    alpha_magic = formulas["alpha_magic"]
    alpha_magic_subs = {
        alpha: alpha_magic,
        diff(alpha, r): diff(alpha_magic, r),
        diff(alpha, r, 2): diff(alpha_magic, r, 2),
    }
    rho_magic = sp.cancel(rho.subs(alpha_magic_subs))
    j_magic = sp.cancel(j_hat.subs(alpha_magic_subs))
    p_r_magic = sp.cancel(p_r.subs(alpha_magic_subs))
    gamma_sq = sp.cancel(1 / (1 - v**2))

    rho_prime_magic = sp.cancel(gamma_sq * (rho_magic - 2 * v * j_magic + v**2 * p_r_magic))
    p_r_prime_magic = sp.cancel(gamma_sq * (v**2 * rho_magic - 2 * v * j_magic + p_r_magic))

    first_integral_magic_rhs = sp.cancel(formulas["first_integral_rhs"].subs(alpha_magic_subs))
    first_integral_magic_expected = sp.cancel(
        8 * pi * (rho - formulas["rho_B"]).subs(alpha_magic_subs) * r**2 * B**3 * alpha_magic
    )

    M_geo_integrand = sp.cancel(8 * pi * formulas["rho_B"] * r**2 * B**3 * alpha_magic)
    M_geo_boundary = sp.cancel(-r**2 * diff(B, r) * (2 + r * diff(B, r) / B))

    checks = [
        check("j_hat{r}(alpha_magic)", j_magic, 0),
        check("p_r(alpha_magic) = -rho(alpha_magic)", p_r_magic, -rho_magic),
        check("rho'(alpha_magic) = rho(alpha_magic)", rho_prime_magic, rho_magic),
        check("p'_r(alpha_magic) = -rho(alpha_magic)", p_r_prime_magic, -rho_magic),
        check("magic-beta first-integral differential form", first_integral_magic_rhs, first_integral_magic_expected),
        check("d/dr M_geo boundary term", sp.cancel(diff(M_geo_boundary, r)), M_geo_integrand),
        check("M_geo = -r^2 B'(1 + alpha_magic)", M_geo_boundary, sp.cancel(-r**2 * diff(B, r) * (1 + alpha_magic))),
        check("det(g) = -alpha^2 det(gamma)", sp.cancel(g.det()), sp.cancel(-alpha**2 * gamma.det())),
        check("g^{tt}", sp.cancel(g_inv[0, 0]), sp.cancel(-1 / alpha**2)),
        check("g^{rr}", sp.cancel(g_inv[1, 1]), sp.cancel(1 / B**2 - beta**2 / alpha**2)),
        check("g^{rr} = -g_tt/(alpha^2 B^2)", sp.cancel(g_inv[1, 1]), sp.cancel(-g[0, 0] / (alpha**2 * B**2))),
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
