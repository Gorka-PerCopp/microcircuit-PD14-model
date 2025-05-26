# -*- coding: utf-8 -*-
#
# network.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-2.0-or-later

#####################
import time
import nest
import numpy as np
import json
from scipy.stats import ks_2samp as ks
import matplotlib.pyplot as plt

## import model implementation
from microcircuit import network
from microcircuit import helpers

## import (default) parameters (network, simulation, stimulus)
from microcircuit.network_params import default_net_dict as net_dict
from microcircuit.sim_params import default_sim_dict as sim_dict
from microcircuit.stimulus_params import default_stim_dict as stim_dict

#####################
populations = net_dict['populations'] # list of populations
t_min = 500.0 # in ms
#####################

## set network scale
scaling_factor = 0.1
net_dict["N_scaling"] = scaling_factor
net_dict["K_scaling"] = scaling_factor

## set path for storing spike data and figures
### TODO revise data path
#sim_dict['data_path'] = '../examples/data_scale_%.2f/' % scaling_factor
seeds = ['12345' + str(i) for i in range(0, 10)] # list of seeds

'''
TODOs:
1. Rate distributions
2. CV distributions
3. CC distributions
'''

#-----------------------------------------------------------------------------------------------------
#                                      RATE DISTRIBUTIONS
#-----------------------------------------------------------------------------------------------------
# load data
rates = {} # list of single neuron firing rates [seed][pop][neuron]
nodes = helpers.json2dict(sim_dict['data_path'] + 'nodes.json')
recording_interval = (max(t_min, sim_dict['t_presim']), sim_dict['t_presim'] + sim_dict['t_sim'])

for cseed, seed in enumerate(seeds):
    data_path = sim_dict['data_path'] + 'seed-%s/' % seed
    rates[cseed] = {}

    for pop in populations:
        rates[cseed][pop] = {}

        label = 'spike_recorder-' + str(nodes['spike_recorder_%s' % pop][0])
        spikes = helpers.load_spike_data(data_path, label)
        rates[cseed][pop] = list(helpers.time_averaged_single_neuron_firing_rates(spikes, nodes[pop], recording_interval))

# store rates as json file

json.dump(rates, open(sim_dict['data_path'] + 'rates.json', 'w'), indent=4)

#------------------------------------------------------------------------------------------------------
#                                      CV DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------
# TODO concatenate CV arrays from all neurons for each population and seed
spike_cvs = {} # list of single neuron CVs [seed][pop][neuron]
for cseed, seed in enumerate(seeds):
    data_path = sim_dict['data_path'] + 'seed-%s/' % seed
    spike_cvs[cseed] = {}

    for pop in populations:
        spike_cvs[cseed][pop] = {}

        label = 'spike_recorder-' + str(nodes['spike_recorder_%s' % pop][0])
        spikes = helpers.load_spike_data(data_path, label)
        spike_cvs[cseed][pop] = list(helpers.single_neuron_isi_cvs(spikes, nodes[pop], recording_interval))

# store spike count variances as json file
json.dump(spike_cvs, open(sim_dict['data_path'] + 'spike_cvs.json', 'w'), indent=4)

#------------------------------------------------------------------------------------------------------
#                                      PAIRWISE-CORRELATION DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------
# load data
cc_binsize = 1. # in ms
spike_ccs = {} # list of pairwise spike count correlations [seed][pop][correlation]
for cseed, seed in enumerate(seeds):
    data_path = sim_dict['data_path'] + 'seed-%s/' % seed
    spike_ccs[cseed] = {}
    
    for pop in populations:
        spike_ccs[cseed][pop] = {}

        label = 'spike_recorder-' + str(nodes['spike_recorder_%s' % pop][0])
        spikes = helpers.load_spike_data(data_path, label)
        spike_ccs[cseed][pop] = list(helpers.pairwise_spike_count_correlations(spikes, nodes[pop], recording_interval, cc_binsize))

# store pairwise spike count correlations as json file
json.dump(spike_ccs, open(sim_dict['data_path'] + 'spike_ccs.json', 'w'), indent=4)

# TODO parameter file for data_analysis
#------------------------------------------------------------------------------------------------------
#                                      ANALYZE RATE DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------
# plot analysis data.py:
# load analysis data
# rates = json.load(open(data_path + 'rates.json', 'r'))

# ## calculate distances (KS) between distributions for different seeds (TODO)
n_seeds = len(seeds)
ks_distances = {} # list of ks distances [pop][seed][ks_distance to other seed]
for cpop, pop in enumerate(populations):
    ks_distances[pop] = {
        "seeds": {},    # values across seeds
        "list": []      # to compute mean and std
    }
    for i in range(n_seeds):
        ks_distances[pop]["seeds"][i] = {}
        for j in range(i+1, n_seeds):
            ks_distance = ks(rates[i][pop], rates[j][pop])[0].tolist()
            #ks_distance = ks(data_hists[i][cpop], data_hists[j][cpop])[0].tolist()
            ks_distances[pop]["seeds"][i][j] = ks_distance
            ks_distances[pop]["list"].append( ks_distance )

# compute the mean KS-distance
mean_ks = []
std_ks = []
for cpop, pop in enumerate(populations):
    ks_list = ks_distances[pop]["list"]#.pop( "list" )
    mean_ks.append(np.mean( ks_list ))
    std_ks.append(np.std( ks_list ))

# save ks distances as json file
json.dump(ks_distances, open(sim_dict['data_path'] + 'ks_distances.json', 'w'), indent=4)
#------------------------------------------------------------------------------------------------------
#                                      PLOTTING
#------------------------------------------------------------------------------------------------------

from matplotlib import rcParams
rcParams['figure.figsize']    = (4,3)
rcParams['figure.dpi']        = 300
rcParams['font.family']       = 'sans-serif'
rcParams['font.size']         = 8
rcParams['legend.fontsize']   = 8
rcParams['axes.titlesize']    = 10
rcParams['axes.labelsize']    = 8
rcParams['ytick.labelsize']   = 8
rcParams['xtick.labelsize']   = 8
rcParams['ytick.major.size']  = 0   ## remove y ticks      
rcParams['text.usetex']       = False 
rcParams['legend.framealpha'] = 1.0
rcParams['legend.edgecolor']  = 'k'

# load analysis data
# ks_distances = json.load(open(data_path + 'ks_distances.json', 'r'))

# compute the best binning for each histogram 
binnings = {} # list of binnings [seed][pop][bins]
for cseed, seed in enumerate(seeds):
    # calculate histogram for each population
    for cpop, pop in enumerate(populations):
        _, bins, _ = helpers.data_distribution(np.array(rates[cseed][pop]), pop, '1/s')
        if cpop not in binnings:
            binnings[cpop] = []
        binnings[cpop].append(bins)

best_bins = {}
for cpop, binning in binnings.items():
    flat_bins = []
    min_diffs = []
    for bins in binning:
        flat_bins += bins.tolist()
        min_diffs.append(np.min(np.diff(bins)))
    max_range = np.max(flat_bins).tolist()
    min_width = np.min(min_diffs).tolist()
    best_bins[cpop] = (max_range, min_width, np.arange(0, max_range + min_width, min_width).tolist())

# calculate histogram for each seed and each population (data_distribution(...))
data_hists = [] # list of histograms [seed][pop][histogram]
data_hist_mat = {}
rate_stats = {} # list of statistics [seed][pop][stats] (mean, std, etc.)
for cseed, seed in enumerate(seeds):
    rate_stats[cseed] = {}
    data_hists.append([])
    for cpop, pop in enumerate(populations):
        rate_stats[cseed][pop] = {}

        rates_pop = np.array(rates[cseed][pop])
        data_hist, bins, stats = helpers.data_distribution(rates_pop, pop, '1/s', np.array(best_bins[cpop][2]))
        data_hists[cseed].append(data_hist.tolist())
        if not cpop in data_hist_mat:
            data_hist_mat[cpop] = np.zeros((len(seeds), len(data_hist)))
        data_hist_mat[cpop][cseed] = data_hist
        rate_stats[cseed][pop] = stats

# save statistics as json file
json.dump(rate_stats, open(sim_dict['data_path'] + 'rate_stats.json', 'w'), indent=4)

# specific example (L23I)
bin_edges_L23I = best_bins[1][2]
bin_centers_L23I = bin_edges_L23I[:-1]  # Left edges for bar alignment
widths_L23I = best_bins[1][1]
bin_err_bar_L23I = [x + widths_L23I / 2 for x in bin_centers_L23I]

means_L23I = np.mean(data_hist_mat[1], axis=0)
stds_L23I = np.std(data_hist_mat[1], axis=0)

plt.figure(1)
plt.bar(bin_centers_L23I, means_L23I, width=widths_L23I, align='edge', edgecolor='black')
plt.errorbar(bin_err_bar_L23I, means_L23I, yerr=stds_L23I, fmt='none', ecolor='red', capsize=5, label='Std Dev')
plt.text(0.35, 0.95, f'mean KS-distance: {mean_ks[1]:.2f}', transform=plt.gca().transAxes, fontsize=12,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

plt.title(r'L23I')
plt.xlabel(r'mean firing rate [$Hz$]')
plt.ylabel(r'frequency')
plt.tight_layout()


plt.figure(2)
grayscale = np.linspace(0.2, 0.8, n_seeds)
colors = [(g, g, g) for g in grayscale]

for cseed in range(n_seeds):
    plt.plot(bin_centers_L23I, data_hist_mat[1][cseed], '-', color=colors[cseed], label=f'Seed {cseed}')

plt.plot(bin_centers_L23I, means_L23I, 'k--', label='Mean')
plt.fill_between(bin_centers_L23I, means_L23I - stds_L23I, means_L23I + stds_L23I, alpha=0.3)
plt.title(r'L23I')
plt.xlabel(r'mean firing rate [$Hz$]')
plt.ylabel(r'frequency')
plt.legend()


# plot of histograms for all populations
fig, axes = plt.subplots(4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.5})

# plot of distributions of KS -distances over seeds
fig1, axes1 = plt.subplots(4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.5})
for cpop, pop in enumerate(populations):
    ax = axes[cpop // 2, cpop % 2]     # axes for each population
    ax1 = axes1[cpop // 2, cpop % 2]

    bin_edges = best_bins[cpop][2]
    bin_centers = bin_edges[:-1]  # Left edges for bar alignment 

    pop_rel_hists = data_hist_mat[cpop] / (scaling_factor * net_dict['full_num_neurons'][cpop])

    # calculate population mean and std histograms across seeds
    pop_mean_hist = np.mean(pop_rel_hists, axis=0)
    pop_std_hist  = np.std(pop_rel_hists, axis=0)

    for cseed in range(n_seeds):
        ax.plot(bin_centers, pop_rel_hists[cseed], '-', color=colors[-1], label=f'Seed {cseed}')

    ax.plot(bin_centers, pop_mean_hist, 'k--', label='Mean')
    ax.fill_between(bin_centers, pop_mean_hist - pop_std_hist, pop_mean_hist + pop_std_hist, alpha=0.3)
    ax.text(0.45, 0.85, f'mean KS-distance: {mean_ks[cpop]:.2f}', transform=ax.transAxes, fontsize=6,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
 
    ax1.hist(ks_distances[pop]["list"], bins=10, color='gray', alpha=0.5)
    ax1.axvline(mean_ks[cpop], color='red', linestyle='--', label='Mean KS-distance')
    ax1.axvline(mean_ks[cpop] + std_ks[cpop], color='blue', linestyle='--', label='Mean + Std')
    ax1.axvline(mean_ks[cpop] - std_ks[cpop], color='blue', linestyle='--', label='Mean - Std')
    ax1.set_xticks(np.arange(0, np.max(ks_distances[pop]["list"]) + 0.1, 0.05))

    ax.set_title(pop)
    ax1.set_title(pop)
    if cpop // 2 == 3:
        ax.set_xlabel(r'mean firing rate [s$^{-1}$]')
        ax1.set_xlabel(r'KS-distance across seeds')
    if cpop % 2 == 0:
        ax.set_ylabel(r'rel. freq.')
        ax1.set_ylabel(r'freq.')
plt.tight_layout()
plt.show()
