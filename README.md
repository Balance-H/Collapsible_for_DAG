# CMDSA Project

This folder contains a cleaned project layout for comparing and benchmarking DAG convex-hull / CMDSA implementations.

## Layout

- `collapsible_nodes_experiment.py` - records average ancestor/collapsed node sizes and their computation time.
- `collapsible_experiment.py` - compares ancestor-model and collapsed-model learning/inference timing and related metrics.
- `data_experiment.py` - complete-data sample-size sweep experiment.
- `csrc/` - C and C++ sources for the `decom_h` Python extension, plus the CMake build file.
- `py/` - standalone Python test scripts and the pure Python DAG convex-hull implementations.
- `scripts/` - helper batch files for building the extension and inspecting dependencies.
- `result/` - output directory for experiment results and CSV files.
- `build/` - generated build artifacts and the compiled `decom_h.pyd` module.

## Experiment Scripts

This project includes three parallel experiment scripts (no primary/secondary distinction):

1. `collapsible_nodes_experiment.py`
	- Focus: structural size and preprocessing cost.
	- Outputs: average `ang_nodes`, average `collapsed_nodes`, and corresponding computation time.
	- Result file: `result/collapsible_nodes_experiment.csv`.

2. `collapsible_experiment.py`
	- Focus: learning/inference performance comparison on selected Bayesian networks.
	- Outputs: aggregated timing/metric comparisons between ancestor-based and collapsed-model pipelines.
	- Result file: `result/collapsible_experiment.csv`.

3. `data_experiment.py`
	- Focus: complete-data setting with sample-size sweep.
	- Design:
	  - No missing-value generation/injection.
	  - Sample sizes: `500, 1000, 2500, 5000`.
	  - Fixed repetitions and settings from script: `r_repeat=50`, `sample_repeat=10`, `r_size=5`.
	  - Network list: `hailfinder`, `hepar2`, `andes`, `munin`.
	  - For each network, the `R` groups are pre-sampled once and reused across all sample sizes to ensure fair comparison.
	- Result file: `result/sample_size_experiment.csv`.

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
