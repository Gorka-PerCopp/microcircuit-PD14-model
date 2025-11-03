# Reference data

## Data generation
[`generate_reference_data.py`](generate_reference_data.py)
Simulate microcircuit with default parameters (in [source](../src/microcircuit/)) for changed simulation time $$T = 10$$s and seed, controlled by [params.py](params.py).
Explain format of simulation data, storage (how and where).

## Data analysis
[analyze_reference_data.py](analyze_reference_data.py)
[compute_ensemble_stats.py](compute_ensemble_stats.py)

## Data visualization
[plot_reference_analysis.py](plot_reference_analysis.py)

## Cluster submission workflow
Specific for slurm queing system.
[`env.sh`](cluster_submission/env.sh)
[`run_benchmark.sh`](cluster_submission/run_benchmark.sh)
[`run_job.sh`](cluster_submission/run_job.sh)


