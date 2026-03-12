# General Relativity Symbolic Calculator

A symbolic General Relativity calculator and report generator powered by SymPy and LaTeX.

## About

This project automates the standard symbolic pipeline in General Relativity so you can move directly from a metric ansatz to geometric and physical quantities without hand-deriving the tensor algebra.

It computes:

- Christoffel symbols
- Riemann curvature tensor
- Ricci tensor and Ricci scalar
- Einstein tensor
- Curvature invariants such as Kretschmann and Weyl
- Orthonormal-frame stress-energy quantities
- Energy conditions
- Geodesic equations
- Killing-coordinate checks
- Bianchi and trace consistency checks

It also generates a LaTeX report and, when `pdflatex` is available, a PDF report.

## Requirements

- Python 3.x
- `sympy`
- a LaTeX distribution with `pdflatex` in PATH

## Quick Start

1. Open [gr_main.py](C:/Users/Nelson/Downloads/GR_python/gr_main.py).
2. In **Section 1**, choose a metric with one line:
   ```python
   METRIC_KEY = 'schwarzschild'
   ```
3. Run:
   ```bash
   python gr_calculator.py
   ```
4. Inspect `gr_report.tex` and `gr_report.pdf`.

## Built-in Metrics

The built-in metric registry currently includes:

- `schwarzschild`
- `reissner_nordstrom`
- `de_sitter_static`
- `minkowski_spherical`
- `frw_flat`
- `static_spherical`
- `morris_thorne_wormhole`
- `pg_areal`
- `pg_spatial_conformal`

These are defined in [gr_metric_library.py](C:/Users/Nelson/Downloads/GR_python/gr_metric_library.py).

## Adding a One-Off Custom Metric

You no longer need to comment and uncomment large metric blocks.

1. In [gr_main.py](C:/Users/Nelson/Downloads/GR_python/gr_main.py), set:
   ```python
   METRIC_KEY = 'custom'
   ```
2. Define any extra parameters or SymPy functions your line element needs in Section 1.2. Typical examples are `Q`, `Lambda`, `Phi(r)`, `b(r)`, or `a(t)`.
3. Fill the `CUSTOM_METRIC_CONFIG` template:
   ```python
   CUSTOM_METRIC_CONFIG = {
       'g_metric': Matrix([...]),
       'metric_name': 'My Custom Metric',
       'metric_description': 'Short description',
       'g_inv_metric': None,
       'e_tetrad': None,
   }
   ```
4. Leave `g_inv_metric` or `e_tetrad` as `None` unless you want to provide them manually.

## Registering a New Built-in Metric

If you want a metric to become a reusable named option:

1. Open [gr_metric_library.py](C:/Users/Nelson/Downloads/GR_python/gr_metric_library.py).
2. Add one more entry inside `build_builtin_metric_library()`.
3. Reuse the same keys as `CUSTOM_METRIC_CONFIG`:
   - `g_metric`
   - `metric_name`
   - `metric_description`
   - optional `g_inv_metric`
   - optional `e_tetrad`
4. Select it later with:
   ```python
   METRIC_KEY = 'your_new_key'
   ```

## Computation Flags

Inside [gr_main.py](C:/Users/Nelson/Downloads/GR_python/gr_main.py), you can still control the expensive steps with flags such as:

- `FAST_MODE = True`
- `COMPUTE_TETRAD = True`
- `COMPUTE_WEYL = True`
- `COMPUTE_KRETSCHMANN = True`
- `OUTPUT_FILENAME = 'gr_report'`

## Project Structure

- [gr_calculator.py](C:/Users/Nelson/Downloads/GR_python/gr_calculator.py): entry point
- [gr_main.py](C:/Users/Nelson/Downloads/GR_python/gr_main.py): user configuration and pipeline
- [gr_metric_library.py](C:/Users/Nelson/Downloads/GR_python/gr_metric_library.py): built-in metric registry
- [gr_tensors.py](C:/Users/Nelson/Downloads/GR_python/gr_tensors.py): symbolic tensor engine
- [gr_latex.py](C:/Users/Nelson/Downloads/GR_python/gr_latex.py): LaTeX/PDF report builder
