# striatal_SPN_lib

This repository contains a library of striatal projection neurons.


How to run the models (would need neuron+python installed, see below)
------------------------------------------------------------------------------

1) compile the mechanisms
2) run example.py


e.g. from a terminal:

    cd mechanisms/single/
    nrnivmodl
    cd ../../
    python3 example.py


Model software
------------------------------------------------------------------------------

NEURON+python: https://www.neuron.yale.edu/neuron/download
(tested on version 7.5 and 7.6)


General
------------------------------------------------------------------------------

Some of the models are used in the publication:

    Lindroos and Hellgren Kotaleski 2020. 
    Predicting complex spikes in striatal projection neurons of the direct pathway 
    following neuromodulation by acetylcholine and dopamine. EJN

Code for simulating, analysing and plotting of ispn is also included (in the ispn branch).
Code for dspn is located in another repo:
https://bitbucket.org/rlindroos/neuron/src/cholinergic_modulation/Complex_spike/


Simulation/analysis/plotting is included in the ispn branch under "Simulations".


Translated versions of the models in: 
    
    Hjorth et al., 2020. 
    The microcircuits of striatum in silico. PNAS
    
are also included. See the included example file and the specific branch for these models.
