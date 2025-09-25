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

default_analysis_dict = {
    # scaling factor of the network
    'scaling_factor': 1.0,
    # RNG seeds
    'RNG_seeds': ['12345' + str(i) for i in range(0, 10)], # RNG_seeds
    # simulation time for analysis
    't_sim': 1.0e+3,
    # start of analysis time interval in ms
    't_min': 500.0,
    # seed for neuron subsampling reroducibility (for CC anlysis)
    'seed_subsampling': 12345, # seed_subsampling
    # subsamples for pairwise statistical analysis
    'subsample_size': 250, # subsample_size
    # bin size for generation of spike-count signals (for CC analysis)
    'binsize': 2.0, # binsize
    # minimal value for CC histogram
    'cc_min': -0.075,
    # maximal value for CC histogram
    'cc_max': 0.1,
    # maximal figure width in inches
    'max_fig_width': 7.5, # max_figure_width 
}
