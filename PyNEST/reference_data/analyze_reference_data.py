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

random.seed( ref_dict['seed_subsampling'] )  # set seed for reproducibility

## set path for storing spike data and figures
### TODO revise data path
#sim_dict['data_path'] = '../examples/data_scale_%.2f/' % scaling_factor
seeds = ref_dict['RNG_seeds'] # list of seeds

sim_dict['data_path'] = 'data_T' + str( int( ref_dict['t_sim'] * 1.0e-3 ) ) + 's/'

########################################################################################################################
#                                   Define auxiliary functions to analyze and plot data                                #
########################################################################################################################

def analyze_single_neuron_stats(observable_name, func):
    observable = {} # list of single neuron observable [seed][pop][neuron]
    recording_interval = ( max( ref_dict['t_min'], ref_dict['t_presim'] ), ref_dict['t_presim'] + ref_dict['t_sim'] )

    for cseed, seed in enumerate( seeds ):
        data_path = sim_dict['data_path'] + 'seed-%s/' % seed
        nodes = helpers.json2dict( data_path + 'nodes.json' )
        observable[cseed] = {}

        for pop in populations:
            observable[cseed][pop] = {}

            label = 'spike_recorder-' + str( nodes['spike_recorder_%s' % pop][0] )
            spikes = helpers.load_spike_data( data_path, label )
            observable[cseed][pop] = list( func( spikes, nodes[pop], recording_interval ) )

    # store observable as json file
    helpers.dict2json( observable, sim_dict['data_path'] + f'{observable_name}.json' )

    return observable

def analyze_pairwise_stats( observable_name, func ):

    recording_interval = ( max ( ref_dict['t_min'], ref_dict['t_presim'] ), ref_dict['t_presim'] + ref_dict['t_sim'] )

    #cc_binsize = 2. # in ms
    observable = {}  # list of pairwise spike count correlations [seed][pop][correlation]
    for cseed, seed in enumerate( seeds ):
        data_path = sim_dict['data_path'] + 'seed-%s/' % seed
        nodes = helpers.json2dict( data_path + 'nodes.json' )
        observable[cseed] = {}
        
        for pop in populations:
            observable[cseed][pop] = {}

            pop_nodes = nodes[pop]  # list of neuron nodes for the population
            label = 'spike_recorder-' + str( nodes['spike_recorder_%s' % pop][0] )
            spikes = helpers.load_spike_data( data_path, label )

            # Get random indices without replacement
            selected_nodes = random.sample( pop_nodes, ref_dict['subsample_size'] )

            observable[cseed][pop] = list( func( spikes, selected_nodes, recording_interval, ref_dict['binsize'] ) )

    # store pairwise spike count correlations as json file
    helpers.dict2json( observable, sim_dict['data_path'] + f'{observable_name}.json' )

    return observable

def compute_ks_distances( observable, observable_name ):
    observable_ks_distances = {} # list of ks distances [pop][seed][ks_distance to other seed]
    for cpop, pop in enumerate( populations ):
        observable_ks_distances[pop] = {
            "seeds": {},    # values across seeds
            "list": []      # to compute mean and std
        }

        # clean data: remove NaNs
        for cseed, seed in enumerate( seeds ):
            observable[cseed][pop] = np.delete( observable[cseed][pop], np.where( np.isnan( observable[cseed][pop] ) ) )

        n_seeds = len( seeds ) 

        for i in range( n_seeds ):
            observable_ks_distances[pop]["seeds"][i] = {}
            for j in range( i+1, n_seeds ):
                observable_ks_distance = ks( observable[i][pop], observable[j][pop] )[0].tolist()
                observable_ks_distances[pop]["seeds"][i][j] = observable_ks_distance
                observable_ks_distances[pop]["list"].append( observable_ks_distance )

    # save ks distances as json file
    helpers.dict2json( observable_ks_distances, sim_dict['data_path'] + f'{observable_name}_ks_distances.json' )

    return observable_ks_distances

def main():
    #cc_binsize = 2.

    rates = analyze_single_neuron_stats( 'rates', helpers.time_averaged_single_neuron_firing_rates )
    spike_cvs = analyze_single_neuron_stats( 'spike_cvs', helpers.single_neuron_isi_cvs )
    spike_ccs = analyze_pairwise_stats( 'spike_ccs', helpers.pairwise_spike_count_correlations )

    compute_ks_distances( rates, 'rate' )
    compute_ks_distances( spike_cvs, 'spike_cvs' )
    compute_ks_distances( spike_ccs, 'spike_ccs' )

    ## current memory consumption of the python process (in MB)
    import psutil
    mem = psutil.Process().memory_info().rss / ( 1024 * 1024 )
    print( f"Current memory consumption: {mem:.2f} MB" )

if __name__ == "__main__":
    main()
