# GR Calculator - Google Colab Edition

Symbolic General Relativity analysis powered by SymPy, with optional GPU numerical post-processing.

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

## Quick Start

1. Open `GR_Colab.ipynb` in Colab.
2. Run Cell 1.
3. Run Cell 2.
4. Edit Cell 3.
5. Run Cell 4 and Cell 5.

## Cell 3 Workflow

Cell 3 is designed to keep metric selection simple.

For the warp-document workflow, the main switches are:

```python
VARIANT = 'variant_a'
PROFILE_MODE = 'document_generic'
RUN_DOCUMENT_COMPARISON = True
GENERATE_COMPARISON_REPORT = True
```

This means:

- `gr_report.pdf` is the direct output of the symbolic run
- `gr_comparison_report.pdf` is optional and only used for checks against external formulas

## Adding a New Metric Outside the Warp Workflow

If you want to work with a brand-new metric beyond the predefined warp variants, use the local driver [gr_main.py](C:/Users/Nelson/Downloads/GR_python/GR_python_colab/gr_main.py) with the same pattern as the desktop version:

1. set `METRIC_KEY = 'custom'`
2. define any extra symbols/functions in Section 1.2
3. fill `CUSTOM_METRIC_CONFIG`

Two commented examples are included directly in [gr_main.py](C:/Users/Nelson/Downloads/GR_python/GR_python_colab/gr_main.py): one diagonal ansatz and one metric with a `dt dr` cross-term.

Reusable built-in metrics are stored in [gr_metric_library.py](C:/Users/Nelson/Downloads/GR_python/GR_python_colab/gr_metric_library.py). Useful examples now include Schwarzschild, Reissner-Nordstrom, static de Sitter, and a Morris-Thorne wormhole template.

## Parallel and GPU Notes

- Symbolic SymPy calculations run on CPU.
- GPU support accelerates only numerical post-processing.
- Parallel symbolic execution is available on Linux/Colab for selected stages.
