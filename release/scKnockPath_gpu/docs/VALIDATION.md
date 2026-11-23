# Release validation

Portability checked on 2026-11-21 on Linux x86_64 with an NVIDIA GPU:

- Extracted the archive outside the repository, then created a **new Conda environment** using `environment.yml` and installed the README's PyTorch 2.8.0 CUDA 12.8 build. Python was 3.11.18.
- Disabled user-site packages and removed inherited Python/dynamic-library paths. NumPy, SciPy, Scanpy, knockpy, groupyr, and PyTorch all imported exclusively from the new environment.
- `python -m pip check` passed with no broken requirements. The installed package versions are listed in [tested-python-packages.txt](tested-python-packages.txt).
- Ran `python tests/smoke_gpu.py` from the extracted directory in that new environment. Both real-data and simulation entry points passed on temporary synthetic inputs (90 cells, 24 genes). Checks independently reconstructed held-out classification errors, the largest-lambda tie rule, W from coefficients, knockoff+ selections, cap flags, and simulation FDP/power.
- Scanned scripts, notebooks, and documentation for machine-specific absolute paths. Runtime model imports resolve relative to the entry point, not the original checkout. Input/output paths are supplied by CLI or are relative to the release root.
- `environment.yml` uses `nodefaults` to avoid inherited default Conda channels and sets `PYTHONNOUSERSITE=1`. The README exports this variable before installation as well, so pre-existing user packages cannot satisfy missing requirements.
- Python syntax/CLI checks and R generator parsing passed. Fig. 2 dry-run listed exactly 510 configurations.
- Model files remain byte-identical to their source copies. Local Markdown links resolve, notebook outputs are empty, and the archive excludes expression data, generated datasets, fitted results, PDFs, and caches.

This validates relocation and a fresh **Python/CUDA** installation on the same Linux host; it is not a physical test on every operating system or GPU. A compatible NVIDIA driver, Conda, network access for dependencies, and sufficient GPU memory remain prerequisites. R dependency installation, full reference downloads, full real-data preprocessing, and manuscript-scale experiments were not rerun during this check.
