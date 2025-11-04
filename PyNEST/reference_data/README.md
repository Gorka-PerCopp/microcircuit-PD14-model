# Reference data

## Data generation
The data generation is exemplified in [`generate_reference_data.py`](generate_reference_data.py).
Here, we simulate the microcircuit with default parameters (in [source](../src/microcircuit/)) for changed simulation time $$T = 10$$s and for ten different RNGseeds, controlled by [`params.py`](params.py).
The data is stored during the simulation process in text files (`.dat`) in the general tree-structure `data/data_T<sim_time_in_s>s/seed-<RNGseed>/`:
* spiking data: `spike_recorder-7717<population-id>-<id>.dat`,
* nodes IDs: `nodes.json`,
* parameters: `sim_dict.json`, `net_dict.json`, and `stim_dict.json`
Explain format of simulation data, storage (how and where).

## Data analysis
[analyze_reference_data.py](analyze_reference_data.py):
* calculates spike statistics (time averaged single neuron firing rates and ISI CVs for each neuron, and spike count correlation coefficients for a subset of 250 neurons) in each population for ten different network realizations.

[compute_ensemble_stats.py](compute_ensemble_stats.py): 
* calculate distributions of spike statistics for each population and each network realization
* calculate Kolmogorov-Smirnoff (KS) distances for each pair of network realizations


## Data visualization
[plot_reference_analysis.py](plot_reference_analysis.py):
* plot distributions of spike statistics for each population and each network realization
* plot distributions of KS distances across pairs of network realizations for each spike statistics and each population

## Cluster submission workflow
Specific for slurm queing system.
[`env.sh`](cluster_submission/env.sh)
[`run_benchmark.sh`](cluster_submission/run_benchmark.sh)
[`run_job.sh`](cluster_submission/run_job.sh)


