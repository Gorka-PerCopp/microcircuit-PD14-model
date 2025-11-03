# Reference data

## Data generation
The data generation is exemplified in [`generate_reference_data.py`](generate_reference_data.py).
Here, we simulate the microcircuit with default parameters (in [source](../src/microcircuit/)) for changed simulation time $$T = 10$$s and for ten different RNGseeds, controlled by [`params.py`](params.py).
The data is stored during the simulation process in text files (`.dat`) in the general tree-structure `data_T<sim_time_in_s>s/seed-<RNGseed>/`:
    * spiking data: `spike_recorder-7717<population-id>-<id>.dat`,
    * rates: `rate<population-id>.dat`
    * nodes IDs: `population_nodeids.dat`
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


