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

## import analysis parameters
from params import params as ref_dict

#####################
populations = net_dict['populations'] # list of populations
#####################

## set network scale
scaling_factor = ref_dict['scaling_factor']
net_dict["N_scaling"] = scaling_factor
net_dict["K_scaling"] = scaling_factor
sim_dict['data_path'] = 'data_T' + str( int( ref_dict['t_sim'] * 1.0e-3 ) ) + 's/'

## set path for storing spike data and figures
### TODO revise data path
#sim_dict['data_path'] = '../examples/data_scale_%.2f/' % scaling_factor
seeds = ref_dict['RNG_seeds'] # list of seeds

########################################################################################################################
#                                   Define auxiliary functions to plot data                                            #
########################################################################################################################

def compute_data_dist( observable, observable_name, units='' ):
    # compute the best binning for each histogram 
    observable_binnings = {} # list of binnings [seed][pop][bins]
    for cseed, seed in enumerate( seeds ):
        cseed = str( cseed )
        # calculate histogram for each population
        for cpop, pop in enumerate( populations ):
            _, bins, _ = helpers.data_distribution( np.array( observable[cseed][pop] ), pop, f'{units}' )
            if cpop not in observable_binnings:
                observable_binnings[cpop] = []
            observable_binnings[cpop].append( bins )

    observable_best_bins = {}
    for cpop, binning in observable_binnings.items():
        max_bins = sorted( binning, key=lambda x: len(x), reverse=True )[0]  # take the binning with the most bins

        max_range = np.max( max_bins ).tolist()
        min_range = np.min( max_bins ).tolist()

        width_diff = np.diff( max_bins )
        min_width = np.min( width_diff ).tolist() if len( width_diff ) > 0 else 0 
        observable_best_bins[cpop] = (min_range, max_range, min_width, np.arange( min_range, max_range + min_width, min_width ).tolist())

    # calculate histogram for each seed and each population (data_distribution(...))
    observable_hists = [] # list of histograms [seed][pop][histogram]
    observable_hist_mat = {}
    observable_stats = {} # list of statistics [seed][pop][stats] (mean, std, etc.)
    for cseed, seed in enumerate( seeds ):
        cseed_str = str( cseed )
        observable_stats[cseed] = {}
        observable_hists.append([])
        for cpop, pop in enumerate( populations ):
            observable_stats[cseed][pop] = {}

            observable_pop = np.array( observable[cseed_str][pop] )
            observable_hist, bins, stats = helpers.data_distribution( observable_pop, pop, f'{units}', np.array( observable_best_bins[cpop][3] ) )
            observable_hists[cseed].append( observable_hist.tolist() )
            if not cpop in observable_hist_mat:
                observable_hist_mat[cpop] = np.zeros( ( len( seeds ), len( observable_hist ) ) )
            observable_hist_mat[cpop][cseed] = observable_hist / stats['sample_size']
            observable_stats[cseed][pop] = stats

    # save statistics as json file
    helpers.dict2json( observable_stats, sim_dict['data_path'] + f'{observable_name}_stats.json' )

    return observable_hist_mat, observable_best_bins, observable_stats

def plot_data_dists(observable_name, x_label, observable_hist_mat, observable_best_bins, observable_ks_distances):

    from matplotlib import rcParams
    rcParams['figure.figsize']    = (ref_dict['max_fig_width'] / 3., (4. / 9.) * ref_dict['max_fig_width'])#(7.5,7)
    rcParams['figure.dpi']        = 300
    rcParams['font.family']       = 'sans-serif'
    rcParams['font.size']         = 8
    rcParams['legend.fontsize']   = 8
    rcParams['axes.titlesize']    = 10
    rcParams['axes.labelsize']    = 8
    rcParams['ytick.labelsize']   = 8
    rcParams['xtick.labelsize']   = 8
    rcParams['ytick.major.size']  = 0   ## remove y ticks      
    rcParams['text.usetex']       = True 
    rcParams['legend.framealpha'] = 1.0
    rcParams['legend.edgecolor']  = 'k'
    data_path = sim_dict['data_path']

    x_max_hist = 0
    x_min_hist = 0
    for cpop, pop in enumerate( populations ):
        bin_edges = observable_best_bins[cpop][3]
        x_max_hist = max( x_max_hist, np.max( bin_edges ) )
        x_min_hist = min( x_min_hist, np.min( bin_edges ) )

    # plot of histograms for all populations
    fig_hist, axes_hist = plt.subplots( 4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.0, 'wspace': 0.0} )

    # plot of distributions of KS -distances
    fig_ks, axes_ks = plt.subplots( 4, 2, sharex=True, sharey=True, gridspec_kw={'hspace': 0.0, 'wspace': 0.0} )

    for cpop, pop in enumerate( populations ):
        ax_hist = axes_hist[cpop // 2, cpop % 2]     # axes for each population
        ax_ks = axes_ks[cpop // 2, cpop % 2]

        bin_edges = observable_best_bins[cpop][3]
        bin_centers = bin_edges[:-1]

        pop_rel_hists = observable_hist_mat[cpop]

        # calculate population mean and std histograms across seeds
        pop_mean_hist = np.mean( pop_rel_hists, axis=0 )
        pop_std_hist  = np.std( pop_rel_hists, axis=0 )

        n_seeds = len( seeds )
        grayscale = np.linspace( 0.2, 0.8, n_seeds )
        colors = [(g, g, g) for g in grayscale]

        for cseed in range( n_seeds ):
            ax_hist.plot( bin_centers, pop_rel_hists[cseed], '-', color=colors[-1], label=f'Seed {cseed}' )

        ax_hist.plot( bin_centers, pop_mean_hist, 'k--', label='Mean' )
        ax_hist.fill_between( bin_centers, pop_mean_hist - pop_std_hist, pop_mean_hist + pop_std_hist, alpha=0.3 )
        
        ks_values = observable_ks_distances[pop]["list"]
        mean = std = None
        if len( ks_values ) > 0:
            mean = np.mean( ks_values )
            std = np.std( ks_values )

        textbox = r'%s' % pop
        if mean is not None:
            textbox += r'\\{\tiny $D_\mathsf{KS} = %.2f$}' % mean
        
            ax_ks.hist( ks_values, bins=n_seeds, color='gray', alpha=0.5 )
            ax_ks.axvline( mean, color='red', linestyle='--', label='Mean KS-distance' )
            ax_ks.axvline( mean + std, color='blue', linestyle='--', label='Mean + Std' )
            ax_ks.axvline( mean - std, color='blue', linestyle='--', label='Mean - Std' )

        ax_hist.text( 0.95, 0.95, textbox, transform=ax_hist.transAxes, fontsize=8,
                verticalalignment='top', horizontalalignment='right' )
        ax_ks.text( 0.95, 0.95, r'%s' % pop, transform=ax_ks.transAxes, fontsize=8,
                verticalalignment='top', horizontalalignment='right' )
        
        ax_ks.set_xlim( 0, np.max( ks_values ) )

        if cpop // 2 == 3:
            ax_ks.set_xlabel( r'KS-distance' )
            
        if cpop % 2 == 0:
            ax_hist.set_ylabel( r'rel. freq.' )
            ax_hist.set_yticks( [] )
            ax_ks.set_ylabel( r'rel. freq.' )
            ax_ks.set_yticks( [] )
        
        if observable_name == 'spike_ccs':
            ax_hist.set_xlim( ref_dict['cc_min'], ref_dict['cc_max'] )
        else:
            ax_hist.set_xlim( 0, x_max_hist )
            ax_hist.set_xticks([0, x_max_hist/2], [r'$0$', r'$%.0f$' % (x_max_hist/2)] )

    fig_hist.text(0.5, -0.01, x_label, ha="center", va="center")

    fig_hist.savefig(f'{data_path}{observable_name}_distributions.pdf',
                 bbox_inches="tight", pad_inches=0.02)
    fig_ks.savefig(f'{data_path}{observable_name}_KS_distances.pdf',
               bbox_inches="tight", pad_inches=0.02)
    

def main():
    data_path = sim_dict['data_path']

    # Read the data in
    rates = helpers.json2dict( f'{data_path}rates.json' )
    spike_cvs = helpers.json2dict( f'{data_path}spike_cvs.json' )
    spike_ccs = helpers.json2dict( f'{data_path}spike_ccs.json' )

    rate_ks_distances = helpers.json2dict( f'{data_path}rate_ks_distances.json' )
    spike_cvs_ks_distances = helpers.json2dict( f'{data_path}spike_cvs_ks_distances.json' )
    spike_ccs_ks_distances = helpers.json2dict( f'{data_path}spike_ccs_ks_distances.json' )

    # Compute distributions and statistics
    rate_hist_mat, rate_best_bins, rate_stats = compute_data_dist( rates, 'rate', '1/s' )
    spike_cvs_hist_mat, spike_cvs_best_bins, spike_cvs_stats = compute_data_dist( spike_cvs, 'spike_cvs' )
    spike_ccs_hist_mat, spike_ccs_best_bins, spike_ccs_stats = compute_data_dist( spike_ccs, 'spike_ccs' )

    # Plot distributions and KS distances
    plot_data_dists( 'rate', r'\begin{center} time averaged single neuron\\firing rate (s$^{-1}$) \end{center}', rate_hist_mat, rate_best_bins, rate_ks_distances )
    plot_data_dists( 'spike_cvs', r'spike irregularity (ISI CV)', spike_cvs_hist_mat, spike_cvs_best_bins, spike_cvs_ks_distances )
    plot_data_dists( 'spike_ccs', r'\begin{center} spike correlation coefficient\\(bin size $%.1f$ ms) \end{center}' % ref_dict['binsize'], spike_ccs_hist_mat, spike_ccs_best_bins, spike_ccs_ks_distances )

    ## current memory consumption of the python process (in MB)
    import psutil
    mem = psutil.Process().memory_info().rss / ( 1024 * 1024 )
    print( f"Current memory consumption: {mem:.2f} MB" )

if __name__ == "__main__":
    main()
    
