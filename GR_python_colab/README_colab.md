# GR Calculator - Google Colab Edition

Symbolic General Relativity analysis powered by SymPy, with optional GPU numerical post-processing.

Spanish quick guide: [GUIA_USO_ES.md](../GUIA_USO_ES.md)

## Repo Intent

This repository is meant to be used in two natural ways:

- as a notebook workflow in Colab
- as a script workflow in Spyder, local Python, or a cloud VM

If you are in Colab, start from `GR_Colab.ipynb`.
If you are local or in Google Cloud, start from `gr_main.py` / `gr_calculator.py` in the project root.

## Files

| File | Purpose |
|---|---|
| `GR_Colab.ipynb` | Main notebook |
| `gr_main.py` | Local driver and configuration |
| `gr_metric_library.py` | Built-in metric registry |
| `gr_tensors.py` | Symbolic tensor engine |
| `gr_latex.py` | LaTeX report assembly |
| `gr_numerics.py` | Numerical evaluation and plotting |
| `gr_warp.py` | Warp-document helpers and comparisons |

## Colab Use Step by Step

1. Open `GR_Colab.ipynb` in Google Colab.
2. Run **Cell 1** to install LaTeX and Python dependencies.
3. Run **Cell 2** to clone the GitHub repo and import the modules.
4. Edit **Cell 3**.
5. Run **Cell 4** for the symbolic calculation.
6. Run **Cell 5** to generate the report files.
7. If you want numerics and plots, continue with **Cell 6** and **Cell 7**. **Cell 8** is an advanced geodesic template and needs a real RHS before it becomes a physical trajectory.

## What Each Important Cell Does

- **Cell 1** installs packages
- **Cell 2** loads the project from GitHub
- **Cell 3** is where you choose the metric/profile/report options
- **Cell 4** runs the symbolic GR computation
- **Cell 5** writes `gr_report.pdf` and, if enabled, `gr_comparison_report.pdf`

## Automatic Tetrad Convention

- `COMPUTE_TETRAD = True` enables the automatic orthonormal frame
- for diagonal metrics, the code uses the canonical coordinate-aligned static tetrad
- for metrics with `g_{0i}` shift terms, it uses the ADM/Eulerian tetrad adapted to the `t = const` slicing
- if you want a different but equivalent tetrad, supply it manually with `e_tetrad = Matrix(...)`

## Standard Warp-Document Workflow

In **Cell 3**, the usual configuration is:

```python
VARIANT = 'variant_a'
PROFILE_MODE = 'document_generic'
RUN_DOCUMENT_COMPARISON = True
GENERATE_COMPARISON_REPORT = True
```

This means:

- `gr_report.pdf` is the direct output of the symbolic run
- `gr_comparison_report.pdf` is optional and only used for checks against external formulas

## If You Want a Brand-New Metric

For a completely new metric outside the predefined warp variants, use the root driver [gr_main.py](../gr_main.py) with the same workflow as the desktop version:

1. set `METRIC_KEY = 'custom'`
2. define any extra symbols/functions in Section 1.2
3. fill `CUSTOM_METRIC_CONFIG`
4. run the script locally, or adapt the same metric definition into a notebook cell

Two commented examples are already included in [gr_main.py](../gr_main.py):

- a diagonal ansatz
- a metric with a `dt dr` cross-term

Reusable built-in metrics are stored in [gr_metric_library.py](../gr_metric_library.py). Useful examples include Schwarzschild, Reissner-Nordstrom, static de Sitter, a Morris-Thorne wormhole template, and the built-in `warp_doc_*` variants from your document.

For a first numerical demo in Colab, switch Cell 3 to `PROFILE_MODE = 'schwarzschild_pg'`.

## GPU and Parallel Notes

- Symbolic SymPy calculations run on CPU.
- GPU support accelerates only numerical post-processing.
- Parallel symbolic execution is available on Linux/Colab for selected stages.
