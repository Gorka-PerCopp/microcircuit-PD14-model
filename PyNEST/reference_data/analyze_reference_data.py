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
import random

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
seed = 12345 # seed for reproducibility
#####################

## set network scale
scaling_factor = 1.
net_dict["N_scaling"] = scaling_factor
net_dict["K_scaling"] = scaling_factor

np.random.seed(seed)  # set seed for reproducibility

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
recording_interval = (max(t_min, sim_dict['t_presim']), sim_dict['t_presim'] + sim_dict['t_sim'])

for cseed, seed in enumerate(seeds):
    data_path = sim_dict['data_path'] + 'seed-%s/' % seed
    nodes = helpers.json2dict(data_path + 'nodes.json')
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
    nodes = helpers.json2dict(data_path + 'nodes.json')
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
spike_ccs = {}  # list of pairwise spike count correlations [seed][pop][correlation]
for cseed, seed in enumerate(seeds):
    data_path = sim_dict['data_path'] + 'seed-%s/' % seed
    nodes = helpers.json2dict(data_path + 'nodes.json')
    spike_ccs[cseed] = {}
    
    for pop in populations:
        spike_ccs[cseed][pop] = {}

        pop_nodes = nodes[pop]  # list of neuron nodes for the population
        label = 'spike_recorder-' + str(nodes['spike_recorder_%s' % pop][0])
        spikes = helpers.load_spike_data(data_path, label)

        # select subset of nodes for the population
        n = len(pop_nodes)
        k = 100

        # Get random indices without replacement
        #selected_indices = random.sample(range(n), k)
        selected_nodes = random.sample(pop_nodes, k)
        '''
        # Get the selected neuron IDs and their spike times
        selected_nodes = [pop_nodes[i] for i in selected_indices]

        # find index of selected nodes in spikes['senders']
        valid_ind = [i for i in selected_nodes if i in spikes['senders']]
        subpop_spikes = {'senders': np.array(valid_ind),
                         'times': spikes['times'][valid_ind]}
        '''
        #spike_ccs[cseed][pop] = list(helpers.pairwise_spike_count_correlations(subpop_spikes, valid_ind, recording_interval, cc_binsize))
        spike_ccs[cseed][pop] = list(helpers.pairwise_spike_count_correlations(spikes, selected_nodes,          recording_interval, cc_binsize))

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
rate_ks_distances = {} # list of ks distances [pop][seed][ks_distance to other seed]
for cpop, pop in enumerate(populations):
    rate_ks_distances[pop] = {
        "seeds": {},    # values across seeds
        "list": []      # to compute mean and std
    }
    for i in range(n_seeds):
        rate_ks_distances[pop]["seeds"][i] = {}
        for j in range(i+1, n_seeds):
            rate_ks_distance = ks(rates[i][pop], rates[j][pop])[0].tolist()
            #ks_distance = ks(data_hists[i][cpop], data_hists[j][cpop])[0].tolist()
            rate_ks_distances[pop]["seeds"][i][j] = rate_ks_distance
            rate_ks_distances[pop]["list"].append( rate_ks_distance )

# compute the mean KS-distance
mean_rate_ks = []
std_rate_ks = []
for cpop, pop in enumerate(populations):
    rate_ks_list = rate_ks_distances[pop]["list"]#.pop( "list" )
    mean_rate_ks.append(np.mean( rate_ks_list ))
    std_rate_ks.append(np.std( rate_ks_list ))

# save ks distances as json file
json.dump(rate_ks_distances, open(sim_dict['data_path'] + 'rate_ks_distances.json', 'w'), indent=4)


#------------------------------------------------------------------------------------------------------
#                                      ANALYZE ISI CV DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------

# spike_cvs = json.load(open(data_path + 'spike_cvs.json', 'r'))
spike_cvs_ks_distances = {} # list of ks distances [pop][seed][ks_distance to other seed]
for cpop, pop in enumerate(populations):
    spike_cvs_ks_distances[pop] = {
        "seeds": {},    # values across seeds
        "list": []      # to compute mean and std
    }
    for i in range(n_seeds):
        spike_cvs_ks_distances[pop]["seeds"][i] = {}
        for j in range(i+1, n_seeds):
            spike_cvs_ks_distance = ks(spike_cvs[i][pop], spike_cvs[j][pop])[0].tolist()
            spike_cvs_ks_distances[pop]["seeds"][i][j] = spike_cvs_ks_distance
            spike_cvs_ks_distances[pop]["list"].append( spike_cvs_ks_distance )

# compute the mean KS-distance
mean_spike_cvs_ks = []
std_spike_cvs_ks = []
for cpop, pop in enumerate(populations):
    spike_cvs_ks_list = spike_cvs_ks_distances[pop]["list"]#.pop( "list" )
    mean_spike_cvs_ks.append(np.mean( spike_cvs_ks_list ))
    std_spike_cvs_ks.append(np.std( spike_cvs_ks_list ))

# save ks distances as json file
json.dump(spike_cvs_ks_distances, open(sim_dict['data_path'] + 'spike_cvs_ks_distances.json', 'w'), indent=4)

#------------------------------------------------------------------------------------------------------
#                                      ANALYZE SPIKE CCS DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------

# spike_cvs = json.load(open(data_path + 'spike_cvs.json', 'r'))
spike_ccs_ks_distances = {} # list of ks distances [pop][seed][ks_distance to other seed]
for cpop, pop in enumerate(populations):
    spike_ccs_ks_distances[pop] = {
        "seeds": {},    # values across seeds
        "list": []      # to compute mean and std
    }
    for i in range(n_seeds):
        spike_ccs_ks_distances[pop]["seeds"][i] = {}
        for j in range(i+1, n_seeds):
            spike_ccs_ks_distance = ks(spike_ccs[i][pop], spike_ccs[j][pop])[0].tolist()
            spike_ccs_ks_distances[pop]["seeds"][i][j] = spike_ccs_ks_distance
            spike_ccs_ks_distances[pop]["list"].append( spike_ccs_ks_distance )

# compute the mean KS-distance
mean_spike_ccs_ks = []
std_spike_ccs_ks = []
for cpop, pop in enumerate(populations):
    spike_ccs_ks_list = spike_ccs_ks_distances[pop]["list"]#.pop( "list" )
    mean_spike_ccs_ks.append(np.mean( spike_ccs_ks_list ))
    std_spike_ccs_ks.append(np.std( spike_ccs_ks_list ))

# save ks distances as json file
json.dump(spike_ccs_ks_distances, open(sim_dict['data_path'] + 'spike_ccs_ks_distances.json', 'w'), indent=4)

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

#------------------------------------------------------------------------------------------------------
#                                      PLOT RATE DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------

# load analysis data
# rate_ks_distances = json.load(open(data_path + 'rate_ks_distances.json', 'r'))

# compute the best binning for each histogram 
rate_binnings = {} # list of binnings [seed][pop][bins]
for cseed, seed in enumerate(seeds):
    # calculate histogram for each population
    for cpop, pop in enumerate(populations):
        _, bins, _ = helpers.data_distribution(np.array(rates[cseed][pop]), pop, '1/s')
        if cpop not in rate_binnings:
            rate_binnings[cpop] = []
        rate_binnings[cpop].append(bins)

rate_best_bins = {}
for cpop, binning in rate_binnings.items():
    flat_bins = []
    min_diffs = []
    for bins in binning:
        flat_bins += bins.tolist()
        min_diffs.append(np.min(np.diff(bins)))
    max_range = np.max(flat_bins).tolist()
    min_width = np.min(min_diffs).tolist()
    rate_best_bins[cpop] = (max_range, min_width, np.arange(0, max_range + min_width, min_width).tolist())

# calculate histogram for each seed and each population (data_distribution(...))
rate_hists = [] # list of histograms [seed][pop][histogram]
rate_hist_mat = {}
rate_stats = {} # list of statistics [seed][pop][stats] (mean, std, etc.)
for cseed, seed in enumerate(seeds):
    rate_stats[cseed] = {}
    rate_hists.append([])
    for cpop, pop in enumerate(populations):
        rate_stats[cseed][pop] = {}

        rates_pop = np.array(rates[cseed][pop])
        rate_hist, bins, stats = helpers.data_distribution(rates_pop, pop, '1/s', np.array(rate_best_bins[cpop][2]))
        rate_hists[cseed].append(rate_hist.tolist())
        if not cpop in rate_hist_mat:
            rate_hist_mat[cpop] = np.zeros((len(seeds), len(rate_hist)))
        rate_hist_mat[cpop][cseed] = rate_hist
        rate_stats[cseed][pop] = stats

# save statistics as json file
json.dump(rate_stats, open(sim_dict['data_path'] + 'rate_stats.json', 'w'), indent=4)

# specific example (L23I)
bin_edges_L23I = rate_best_bins[1][2]
bin_centers_L23I = bin_edges_L23I[:-1]  # Left edges for bar alignment
widths_L23I = rate_best_bins[1][1]
bin_err_bar_L23I = [x + widths_L23I / 2 for x in bin_centers_L23I]

means_L23I = np.mean(rate_hist_mat[1], axis=0)
stds_L23I = np.std(rate_hist_mat[1], axis=0)

plt.figure(1)
plt.bar(bin_centers_L23I, means_L23I, width=widths_L23I, align='edge', edgecolor='black')
plt.errorbar(bin_err_bar_L23I, means_L23I, yerr=stds_L23I, fmt='none', ecolor='red', capsize=5, label='Std Dev')
plt.text(0.35, 0.95, f'mean KS-distance: {mean_rate_ks[1]:.2f}', transform=plt.gca().transAxes, fontsize=12,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

plt.title(r'L23I')
plt.xlabel(r'mean firing rate [$Hz$]')
plt.ylabel(r'frequency')
plt.tight_layout()


plt.figure(2)
grayscale = np.linspace(0.2, 0.8, n_seeds)
colors = [(g, g, g) for g in grayscale]

for cseed in range(n_seeds):
    plt.plot(bin_centers_L23I, rate_hist_mat[1][cseed], '-', color=colors[cseed], label=f'Seed {cseed}')

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

    bin_edges = rate_best_bins[cpop][2]
    bin_centers = bin_edges[:-1]  # Left edges for bar alignment 

    pop_rel_hists = rate_hist_mat[cpop] / (scaling_factor * net_dict['full_num_neurons'][cpop])

    # calculate population mean and std histograms across seeds
    pop_mean_hist = np.mean(pop_rel_hists, axis=0)
    pop_std_hist  = np.std(pop_rel_hists, axis=0)

    for cseed in range(n_seeds):
        ax.plot(bin_centers, pop_rel_hists[cseed], '-', color=colors[-1], label=f'Seed {cseed}')

    ax.plot(bin_centers, pop_mean_hist, 'k--', label='Mean')
    ax.fill_between(bin_centers, pop_mean_hist - pop_std_hist, pop_mean_hist + pop_std_hist, alpha=0.3)
    ax.text(0.45, 0.85, f'mean KS-distance: {mean_rate_ks[cpop]:.2f}', transform=ax.transAxes, fontsize=6,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
 
    ax1.hist(rate_ks_distances[pop]["list"], bins=10, color='gray', alpha=0.5)
    ax1.axvline(mean_rate_ks[cpop], color='red', linestyle='--', label='Mean KS-distance')
    ax1.axvline(mean_rate_ks[cpop] + std_rate_ks[cpop], color='blue', linestyle='--', label='Mean + Std')
    ax1.axvline(mean_rate_ks[cpop] - std_rate_ks[cpop], color='blue', linestyle='--', label='Mean - Std')
    ax1.set_xticks(np.arange(0, np.max(rate_ks_distances[pop]["list"]) + 0.1, 0.05))

    ax.set_title(pop)
    ax1.set_title(pop)
    if cpop // 2 == 3:
        ax.set_xlabel(r'mean firing rate [s$^{-1}$]')
        ax1.set_xlabel(r'KS-distance across seeds')
    if cpop % 2 == 0:
        ax.set_ylabel(r'rel. freq.')
        ax1.set_ylabel(r'freq.')
plt.tight_layout()
plt.savefig('data/rate_distributions.pdf')

#------------------------------------------------------------------------------------------------------
#                                      PLOT ISI CV DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------

# load analysis data
# spike_cvs_ks_distances = json.load(open(data_path + 'spike_cvs_ks_distances.json', 'r'))

# compute the best binning for each histogram 
spike_cvs_binnings = {} # list of binnings [seed][pop][bins]
for cseed, seed in enumerate(seeds):
    # calculate histogram for each population
    for cpop, pop in enumerate(populations):
        _, bins, _ = helpers.data_distribution(np.array(spike_cvs[cseed][pop]), pop, '')
        if cpop not in spike_cvs_binnings:
            spike_cvs_binnings[cpop] = []
        spike_cvs_binnings[cpop].append(bins)

spike_cvs_best_bins = {}
for cpop, binning in spike_cvs_binnings.items():
    flat_bins = []
    min_diffs = []
    for bins in binning:
        flat_bins += bins.tolist()
        min_diffs.append(np.min(np.diff(bins)))
    max_range = np.max(flat_bins).tolist()
    min_width = np.min(min_diffs).tolist()
    spike_cvs_best_bins[cpop] = (max_range, min_width, np.arange(0, max_range + min_width, min_width).tolist())

# calculate histogram for each seed and each population (data_distribution(...))
spike_cvs_hists = [] # list of histograms [seed][pop][histogram]
spike_cvs_hist_mat = {}
spike_cvs_stats = {} # list of statistics [seed][pop][stats] (mean, std, etc.)
for cseed, seed in enumerate(seeds):
    spike_cvs_stats[cseed] = {}
    spike_cvs_hists.append([])

    for cpop, pop in enumerate(populations):
        spike_cvs_stats[cseed][pop] = {}

        spike_cvs_pop = np.array(spike_cvs[cseed][pop])
        spike_cvs_hist, bins, stats = helpers.data_distribution(spike_cvs_pop, pop, '', np.array(spike_cvs_best_bins[cpop][2]))
        spike_cvs_hists[cseed].append(spike_cvs_hist.tolist())

        if not cpop in spike_cvs_hist_mat:
            spike_cvs_hist_mat[cpop] = np.zeros((len(seeds), len(spike_cvs_hist)))
        spike_cvs_hist_mat[cpop][cseed] = spike_cvs_hist / stats['sample_size']
        spike_cvs_stats[cseed][pop] = stats

# plot of histograms for all populations
fig2, axes2 = plt.subplots(4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.5})

# plot of distributions of KS -distances over seeds
fig3, axes3 = plt.subplots(4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.5})
for cpop, pop in enumerate(populations):
    ax2 = axes2[cpop // 2, cpop % 2]     # axes for each population
    ax3 = axes3[cpop // 2, cpop % 2]

    bin_edges = spike_cvs_best_bins[cpop][2]
    bin_centers = bin_edges[:-1]  # Left edges for bar alignment 

    #for cseed in range(n_seeds):
    pop_rel_hists = spike_cvs_hist_mat[cpop]

    # calculate population mean and std histograms across seeds
    pop_mean_hist = np.mean(pop_rel_hists, axis=0)
    pop_std_hist  = np.std(pop_rel_hists, axis=0)

    for cseed in range(n_seeds):
        ax2.plot(bin_centers, pop_rel_hists[cseed], '-', color=colors[-1], label=f'Seed {cseed}')

    ax2.plot(bin_centers, pop_mean_hist, 'k--', label='Mean')
    ax2.fill_between(bin_centers, pop_mean_hist - pop_std_hist, pop_mean_hist + pop_std_hist, alpha=0.3)
    ax2.text(0.45, 0.85, f'mean KS-distance: {mean_spike_cvs_ks[cpop]:.2f}', transform=ax2.transAxes, fontsize=6,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
 
    ax3.hist(spike_cvs_ks_distances[pop]["list"], bins=10, color='gray', alpha=0.5)
    ax3.axvline(mean_spike_cvs_ks[cpop], color='red', linestyle='--', label='Mean KS-distance')
    ax3.axvline(mean_spike_cvs_ks[cpop] + std_spike_cvs_ks[cpop], color='blue', linestyle='--', label='Mean + Std')
    ax3.axvline(mean_spike_cvs_ks[cpop] - std_spike_cvs_ks[cpop], color='blue', linestyle='--', label='Mean - Std')
    ax3.set_xticks(np.arange(0, np.max(spike_cvs_ks_distances[pop]["list"]) + 0.1, 0.05))

    ax2.set_title(pop)
    ax3.set_title(pop)
    if cpop // 2 == 3:
        ax2.set_xlabel(r'ISI CV')
        ax3.set_xlabel(r'KS-distance across seeds')
    if cpop % 2 == 0:
        ax2.set_ylabel(r'rel. freq.')
        ax3.set_ylabel(r'freq.')
plt.tight_layout()
plt.savefig('data/ISI_CV_distributions.pdf')

#------------------------------------------------------------------------------------------------------
#                                      PLOT SPIKE CCS DISTRIBUTIONS
#------------------------------------------------------------------------------------------------------

# load analysis data
# spike_ccs_ks_distances = json.load(open(data_path + 'spike_ccs_ks_distances.json', 'r'))

# compute the best binning for each histogram 
spike_ccs_binnings = {} # list of binnings [seed][pop][bins]
for cseed, seed in enumerate(seeds):
    # calculate histogram for each population
    for cpop, pop in enumerate(populations):
        _, bins, _ = helpers.data_distribution(np.array(spike_ccs[cseed][pop]), pop, '')
        if cpop not in spike_ccs_binnings:
            spike_ccs_binnings[cpop] = []
        spike_ccs_binnings[cpop].append(bins)


spike_ccs_best_bins = {}
for cpop, binning in spike_ccs_binnings.items():
    max_bins = sorted(binning, key=lambda x: len(x), reverse=True)[0]  # take the binning with the most bins
    max_range = np.max(max_bins).tolist()
    width_diff = np.diff(max_bins)
    min_width = np.min(width_diff).tolist() if len(width_diff) > 0 else 0  # avoid division by zero
    spike_ccs_best_bins[cpop] = (max_range, min_width, max_bins.tolist())
    print(f"Expected bin length for pop {cpop}: {len(spike_ccs_best_bins[cpop][2])}")

# calculate histogram for each seed and each population (data_distribution(...))
spike_ccs_hists = [] # list of histograms [seed][pop][histogram]
spike_ccs_hist_mat = {}
spike_ccs_stats = {} # list of statistics [seed][pop][stats] (mean, std, etc.)
for cseed, seed in enumerate(seeds):
    spike_ccs_stats[cseed] = {}
    spike_ccs_hists.append([])

    for cpop, pop in enumerate(populations):
        spike_ccs_stats[cseed][pop] = {}

        spike_ccs_pop = np.array(spike_cvs[cseed][pop])
        spike_ccs_hist, bins, stats = helpers.data_distribution(spike_ccs_pop, pop, '', np.array(spike_ccs_best_bins[cpop][2]))
        spike_ccs_hists[cseed].append(spike_ccs_hist.tolist())

        if not cpop in spike_ccs_hist_mat:
            spike_ccs_hist_mat[cpop] = np.zeros((len(seeds), len(spike_ccs_hist)))
        spike_ccs_hist_mat[cpop][cseed] = spike_ccs_hist / stats['sample_size']
        spike_ccs_stats[cseed][pop] = stats
'''
# plot of histograms for all populations
fig3, axes3 = plt.subplots(4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.5})

# plot of distributions of KS -distances over seeds
fig4, axes4 = plt.subplots(4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.5})
for cpop, pop in enumerate(populations):
    ax3 = axes3[cpop // 2, cpop % 2]     # axes for each population
    ax4 = axes4[cpop // 2, cpop % 2]

    bin_edges = spike_ccs_best_bins[cpop][2]
    bin_centers = bin_edges[:-1]  # Left edges for bar alignment 

    #for cseed in range(n_seeds):
    pop_rel_hists = spike_ccs_hist_mat[cpop]

    # calculate population mean and std histograms across seeds
    pop_mean_hist = np.mean(pop_rel_hists, axis=0)
    pop_std_hist  = np.std(pop_rel_hists, axis=0)

    for cseed in range(n_seeds):
        ax2.plot(bin_centers, pop_rel_hists[cseed], '-', color=colors[-1], label=f'Seed {cseed}')

    ax2.plot(bin_centers, pop_mean_hist, 'k--', label='Mean')
    ax2.fill_between(bin_centers, pop_mean_hist - pop_std_hist, pop_mean_hist + pop_std_hist, alpha=0.3)
    ax2.text(0.45, 0.85, f'mean KS-distance: {mean_spike_ccs_ks[cpop]:.2f}', transform=ax2.transAxes, fontsize=6,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
 
    ax3.hist(spike_ccs_ks_distances[pop]["list"], bins=10, color='gray', alpha=0.5)
    ax3.axvline(mean_spike_ccs_ks[cpop], color='red', linestyle='--', label='Mean KS-distance')
    ax3.axvline(mean_spike_ccs_ks[cpop] + std_spike_ccs_ks[cpop], color='blue', linestyle='--', label='Mean + Std')
    ax3.axvline(mean_spike_ccs_ks[cpop] - std_spike_ccs_ks[cpop], color='blue', linestyle='--', label='Mean - Std')
    ax3.set_xticks(np.arange(0, np.max(spike_ccs_ks_distances[pop]["list"]) + 0.1, 0.05))

    ax2.set_title(pop)
    ax3.set_title(pop)
    if cpop // 2 == 3:
        ax2.set_xlabel(r'CC')
        ax3.set_xlabel(r'KS-distance across seeds')
    if cpop % 2 == 0:
        ax2.set_ylabel(r'rel. freq.')
        ax3.set_ylabel(r'freq.')
plt.tight_layout()
'''
