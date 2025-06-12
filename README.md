# Cortical microcircuit model (Potjans & Diesmann, 2014)

## Overview

This [repository](https://github.com/INM-6/microcircuit-PD14-model) contains a detailed mathematical description and a reference implementation of the model of a cortical microcircuit proposed by [Potjans & Diesmann (2014)][1]. 
The model describes the neuronal circuitry under ~1 mm$`^2`$ cortical surface. 
It comprises four cortical layers (L2/3, L4, L5, L6), each represented by a randomly connected network of excitatory and inhibitory point neurons. 
The network connectivity is derived from anatomical and electrophysiological data.
Connection probabilities between neurons in the network are highly specific and depend on the cell type (excitatory, inhibitory) and on the locations (cortical layers) of the pre- and postsynaptic neurons.
In contrast to this high specificity in the connectivity, all neurons in the network are identical and share the same dynamics and parameters, irrespective of their type and location.
Similarly, all synapses are described by an identical dynamics, and differ only in the synaptic weight and spike-transmission latencies.
Synaptic weights and spike transmission latencies are randomly drawn from distributions which depend only on the type of the presynaptic neuron (excitatory or inhibitory), but are otherwise identical for all neurons and connections (with one exception).
In addition to inputs from the local network, neurons receive external inputs representing thalamic afferents and cortico-cortical inputs from more distant cortical regions. 

The original purpose of this model was to understand the relationship between the connectivity and the spiking activity within local cortical circuits. 
Specifically, the model demonstrates that the observed cell-type and layer specificity of in-vivo firing rates is largely explained by the specificity in the number of connections between cortical subpopulations, and doesn't require a specificity in single neuron or synapse dynamics.

|  |  |  |
|--|--|--|
| <img src="figures/potjans_2014_microcircuit.png" width="300"/> | <img src="figures/potjans_2014_raster_plot.png" width="400"/> | <img src="figures/potjans_2014_box_plot.png" width="400"/> |

*Sketch of the cortical microcircuit model (left), spiking activity (middle) and distributions of time averaged single-neuron firing rates across neurons in each subpopulation (right). Adapted from ([van Albada et al., 2018][2])*

In recent years, the model became an established Computational Neuroscience benchmark for various soft- and hardware architectures ([van Albada et al., 2018][2]; [Jordan et al., 2018][3]; [Rhodes et al., 2020][4]; [Dasbach et al., 2021][5]; [Albers et al., 2022][6]; [Kurth et al., 2022][7]; [Heittmann et al., 2022][8]; [Pronold et al., 2022][9]; [Pronold et al., 2022][10]; [Golosio et al., 2023][11]; [Kauth et al., 2023][12]; [Schmidt et al., 2024][13]). 

A list of studies citing and/or using the microcircuit model is provided [here](publications.bib) (studies *using* the model are tagged by `uses_PD14 = {yes}`). Feel free to contact us in case publications are missing in this list.

## Model description

A detailed mathematical, implementation agnostic description of the model and its parameters is provided [here](docs/ModelDescription_microcircuit-PD14-model.pdf).

## Model implementations
* [PyNEST](PyNEST/README.md)

## Repository contents

|  |  | 
|--|--|
| [`docs`](docs) | model description (implementation agnostic)|
| [`PyNEST`](PyNEST) | PyNEST implementation (python package)|
| &emsp;[`PyNEST/src/microcircuit`](PyNEST/src/microcircuit) | source code |
| &emsp;[`PyNEST/examples`](PyNEST/examples) | examples illustrating usage of the python package |
| &emsp;[`PyNEST/reference_data`](PyNEST/reference_data) | reference spike data |
| &emsp;[`PyNEST/tests`](PyNEST/tests) | unit tests |
| [`publications.bib`](publications.bib) | publications citing/using the microcircuit model |
| [`figures`](figures) | overview figures |

## References

[1]: <https://doi.org/10.1093/cercor/bhs358> "Potjans & Diesmann (2014). The cell-type specific cortical microcircuit: relating structure and activity in a full-scale spiking network model. Cerebral Cortex, 24(3), 785-806."
[Potjans & Diesmann (2014). The cell-type specific cortical microcircuit: relating structure and activity in a full-scale spiking network model. Cerebral Cortex, 24(3), 785-806.](https://doi.org/10.1093/cercor/bhs358)

[2]: <https://doi.org/10.3389/fnins.2018.00291> "van Albada et al. (2018). Performance Comparison of the Digital Neuromorphic Hardware SpiNNaker and the Neural Network Simulation Software NEST for a Full-Scale Cortical Microcircuit Model. Front. Neurosci. 12:291."
[van Albada et al. (2018). Performance Comparison of the Digital Neuromorphic Hardware SpiNNaker and the Neural Network Simulation Software NEST for a Full-Scale Cortical Microcircuit Model. Front. Neurosci. 12:291.](https://doi.org/10.3389/fnins.2018.00291)

[3]: <https://doi.org/10.3389/fninf.2018.00002> "Jordan et al. (2018). Extremely scalable spiking neuronal network simulation code: from laptops to exascale computers. Frontiers in Neuroinformatics, 12, 317068."
[Jordan et al. (2018). Extremely scalable spiking neuronal network simulation code: from laptops to exascale computers. Frontiers in Neuroinformatics, 12, 317068.](https://doi.org/10.3389/fninf.2018.00002)

[4]: <https://doi.org/10.1098/rsta.2019.0160> "Rhodes et al. (2020). Real-time cortical simulation on neuromorphic hardware. Philosophical Transactions of the Royal Society A, 378(2164), 20190160."
[Rhodes et al. (2020). Real-time cortical simulation on neuromorphic hardware. Philosophical Transactions of the Royal Society A, 378(2164), 20190160.](https://doi.org/10.1098/rsta.2019.0160)

[5]: <https://doi.org/10.3389/fnins.2021.757790> "Dasbach et al. (2021). Dynamical characteristics of recurrent neuronal networks are robust against low synaptic weight resolution. Frontiers in Neuroscience, 15, 757790."
[Dasbach et al. (2021). Dynamical characteristics of recurrent neuronal networks are robust against low synaptic weight resolution. Frontiers in Neuroscience, 15, 757790.](https://doi.org/10.3389/fnins.2021.757790)

[6]: <https://doi.org/10.3389/fninf.2022.837549> "Albers et al. (2022). A modular workflow for performance benchmarking of neuronal network simulations. Frontiers in neuroinformatics, 16, 837549."
[Albers et al. (2022). A modular workflow for performance benchmarking of neuronal network simulations. Frontiers in neuroinformatics, 16, 837549.](https://doi.org/10.3389/fninf.2022.837549)

[7]: <https://doi.org/10.1088/2634-4386/ac55fc> "Kurth et al. (2022). Sub-realtime simulation of a neuronal network of natural density. Neuromorphic computing and engineering, 2(2), 021001."
[Kurth et al. (2022). Sub-realtime simulation of a neuronal network of natural density. Neuromorphic computing and engineering, 2(2), 021001.](https://doi.org/10.1088/2634-4386/ac55fc)

[8]: <https://doi.org/10.3389/fnins.2021.728460> "Heittmann et al. (2022). Simulating the cortical microcircuit significantly faster than real time on the IBM INC-3000 neural supercomputer. Frontiers in Neuroscience, 15, 728460."
[Heittmann et al. (2022). Simulating the cortical microcircuit significantly faster than real time on the IBM INC-3000 neural supercomputer. Frontiers in Neuroscience, 15, 728460.](https://doi.org/10.3389/fnins.2021.728460)

[9]: <https://doi.org/10.3389/fninf.2021.785068> "Pronold et al. (2022). Routing brain traffic through the von Neumann bottleneck: Parallel sorting and refactoring. Frontiers in Neuroinformatics, 15, 785068."
[Pronold et al. (2022). Routing brain traffic through the von Neumann bottleneck: Parallel sorting and refactoring. Frontiers in Neuroinformatics, 15, 785068.](https://doi.org/10.3389/fninf.2021.785068)

[10]: <https://doi.org/10.1016/j.parco.2022.102952> "Pronold et al. (2022). Routing brain traffic through the von Neumann bottleneck: Efficient cache usage in spiking neural network simulation code on general purpose computers. Parallel computing, 113, 102952."
[Pronold et al. (2022). Routing brain traffic through the von Neumann bottleneck: Efficient cache usage in spiking neural network simulation code on general purpose computers. Parallel computing, 113, 102952.](https://doi.org/10.1016/j.parco.2022.102952)

[11]: <https://doi.org/10.3390/app13179598> "Golosio et al.  (2023). Runtime Construction of Large-Scale Spiking Neuronal Network Models on GPU Devices. Applied Sciences, 13(17), 9598."
[Golosio et al.  (2023). Runtime Construction of Large-Scale Spiking Neuronal Network Models on GPU Devices. Applied Sciences, 13(17), 9598.](https://doi.org/10.3390/app13179598)

[12]: <https://doi.org/10.3389/fncom.2023.1144143> "Kauth et al. (2023). neuroAIx-Framework: design of future neuroscience simulation systems exhibiting execution of the cortical microcircuit model 20× faster than biological real-time. Frontiers in Computational Neuroscience, 17, 1144143."
[Kauth et al. (2023). neuroAIx-Framework: design of future neuroscience simulation systems exhibiting execution of the cortical microcircuit model 20× faster than biological real-time. Frontiers in Computational Neuroscience, 17, 1144143.](https://doi.org/10.3389/fncom.2023.1144143)

[13]: <https://doi.org/10.48550/arXiv.2412.02619> "Schmidt et al. (2024). Demonstrating the Advantages of Analog Wafer-Scale Neuromorphic Hardware. arXiv:2412.02619."
[Schmidt et al. (2024). Demonstrating the Advantages of Analog Wafer-Scale Neuromorphic Hardware. arXiv:2412.02619.](https://doi.org/10.48550/arXiv.2412.02619)

## Contact
- [Johanna Senk](mailto:j.senk@fz-juelich.de)
- [Tom Tetzlaff](mailto:t.tetzlaff@fz-juelich.de)

## Contribute
We welcome contributions to the documentation and the code. For bug reports, feature requests, documentation improvements, or other issues, please create a [GitHub issue](https://github.com/INM-6/microcircuit-PD14-model/issues/new/choose).

## License

The material in this repository is subject to different licenses:

- All material outside the `PyNEST` folder is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa]. For details, see [here](LICENSES/CC-BY-NC-SA-4.0.txt). 
  [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

- The material inside the `PyNEST` folder is licensed under the [GNU General Public License v2.0 or later](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html). For details, see [here](LICENSES/GPL-2.0-or-later.txt). 
  [![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

