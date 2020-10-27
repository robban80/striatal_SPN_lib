# striatal_SPN_lib

This is the ispn branch of the repo.

This readme is UNDER CONSTRUCTION...

The example file included here does not add anything to the ones in the main branch 
(and will probably be removed)
By default it simulates the five "first" models of the ispn lib - giving a current injection of random amplitude (rheobase + 0-60 pA).

The files used for simulation and analysis can be found under "Simulations".


Simulations
-----------------------------------------------------------------------------
Simulations were run the Beskow super computer, PDC-KTH:
https://www.pdc.kth.se/hpc-services/computing-systems/beskow-1.737436

* ispn_run_inVivo_ramping_randMod.py was used to run the simulations
* the file.job was used to start the simulations
* functions4analysis.py and plot.py were used to plot subpanels for fig 6 in Lindroos & H. K. (2020) (there is a similar file in the dspn repo---se "Genereal", below)

Plotting subpanels
------------------------------------------------------------------------------

The file Simulations/plot.py can be used to replot the panels of figure 6 in Lindroos & H.K. (2020).
This is done using aggregated data (in .json files). The raw data will be uploaded to OSF:
https://osf.io/wmsdj/

To replot the data, do the following in a terminal:
cd Simulations
python3 plot.py

How to run the models (would need neuron+python installed, see below)
------------------------------------------------------------------------------
To run the simulation script locally:

cd Simulatins
nrnivmodl ../mechanisms/single
python3 ispn_run_inVivo_ramping_randMod.py


Model software
------------------------------------------------------------------------------

NEURON+python: https://www.neuron.yale.edu/neuron/download
(tested on version 7.5 and 7.6)

python3
mpi4py (to run in parallel)
sys, numpy, json, pickle
(some additional packages are requested by plot.py. Of these matplotlib
is definitely needed, but the others can likely be commented)

General (other models)
------------------------------------------------------------------------------

Code for dspn is located in another repo:
https://bitbucket.org/rlindroos/neuron/src/cholinergic_modulation/Complex_spike/
(in branch: cholinergic_modulation)


Translated versions of the models in: 
    
    Hjorth et al., 2020. 
    The microcircuits of striatum in silico. PNAS
    
are also included (in the simulating_network_models branch ...and main?). See the included example file and the specific branch for these models.

