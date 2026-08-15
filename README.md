# CMDSA Project

This folder contains a cleaned project layout for comparing and benchmarking DAG convex-hull / CMDSA implementations.

## Layout

- `collapsible_experiment.py` - final end-to-end experiment script.
- `collapsible_nodes_experiment.py` - lightweight timing script that records average ancestor/collapsed node sizes and computation time.
- `csrc/` - C and C++ sources for the `decom_h` Python extension, plus the CMake build file.
- `py/` - standalone Python test scripts and the pure Python DAG convex-hull implementations.
- `scripts/` - helper batch files for building the extension and inspecting dependencies.
- `result/` - output directory for experiment results and CSV files.
- `build/` - generated build artifacts and the compiled `decom_h.pyd` module.

## Main Experiment

Run `collapsible_experiment.py` to compare ancestor-model and collapsed-model learning/inference times on the selected Bayesian network examples. The script writes aggregated results into `result/collapsible_experiment.csv`.

Run `collapsible_nodes_experiment.py` when you only need the average `ang_nodes` size, average `collapsed_nodes` size, and the time spent computing each of them. The script writes aggregated results into `result/collapsible_nodes_experiment.csv`.

## Dependencies

### Python

- Python 3.9
- pandas
- networkx
- numpy
- igraph
- pgmpy

### Native

- Microsoft Visual Studio C++ build tools
- CMake
- pybind11
- igraph native libraries
- uthash header

## Build

Use `scripts/build_pyd.cmd` for a one-click rebuild of the `decom_h.pyd` extension.
`scripts/build_decom_h.cmd` is still available if you prefer the more explicit rebuild flow that clears the cache first.

## Notes

- The Python scripts expect the compiled extension to be in `build/`.
- `result/` is intended for experiment outputs only.
- The project is organized so the root folder keeps only the final experiment entry point and project directories.