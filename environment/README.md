# Reproducible Python environment

Python 3.11 is the manuscript baseline. Python 3.12 is supported only after
both Python versions pass the same unit, integration, and scientific snapshot
tests.

Create and activate the public conda-forge environment from the repository
root:

```console
conda env create --file environment/environment.yml
conda activate oips-repro
python -m pip install --editable .
```
