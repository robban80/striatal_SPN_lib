
# script used for running in vivo simultion of neuromod
# This script is based on plateaus_in_vivo_saveFeatures.py -> inVivo_run_mixedPattern.py
#
# Results are stored under Results/InVivo_sim_res/InVivo_clustered_D<>_section<>.pkl

from neuron import h

import sys
sys.path.insert(0, '../')

import json, pickle
import numpy                as np
import MSN_builder          as build 
import common_functions     as use
import modulation_lib       as modulate

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

pc = h.ParallelContext()
id = int(pc.id())

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# ----------------------------------------------------------------------------------------

def run_model(cell_index, ci):
    
    # 1s low input; 0.3s ramping input; 0.2s high input
    tstop       = 1500
    modOnTime   = 1000
    nbg         = 5
    tau         = 300
    conditions  = ['ACh+DA','ACh','DA','ctrl']
    c           = conditions[ci]
    
    # channel distribution parameter library   
    path            = '../Libraries'
    with open('{}/D2_34bestFit_updRheob.pkl'.format(path), 'rb') as f:
        model_sets  = pickle.load(f, encoding="latin1")

    parameters      =   model_sets[cell_index]['variables'] 
    par             =   '../params_iMSN.json'
    morphology      =   '../Morphologies/WT-dMSN_P270-20_1.02_SGA1-m24.swc'

    cell = build.MSN(  params=par,                  \
                       morphology=morphology,       \
                       variables=parameters         )

    # record vectors
    tm  = h.Vector()
    vm  = h.Vector()
    tm.record(h._ref_t)
    vm.record(cell.soma(0.5)._ref_v)
    
    # transient to play
    transient = h.Vector([use.alpha(ht, modOnTime, 1, tau) if ht >= modOnTime else 0 for ht in np.arange(0,tstop,h.dt)])
    
    
    # result dict
    res = {}
    
    for i in range(20):
        res[i] = {}
        # draw random factors
        modulation_DA = {'intr':   {'naf': np.random.uniform(0.95,1.1),
                                    'kaf': np.random.uniform(1.0,1.1),
                                    'kas': np.random.uniform(1.0,1.1),
                                    'kir': np.random.uniform(0.8,1.0),
                                    'can': np.random.uniform(0.9,1.0),
                                    'car': np.random.uniform(0.6,0.8),
                                    'cal12': np.random.uniform(0.7,0.8),
                                    'cal13': np.random.uniform(0.7,0.8)
                                    },
                         'syn':    {'NMDA':np.random.uniform(0.85,1.05), 
                                    'AMPA':np.random.uniform(0.7,0.9), 
                                    'GABA':np.random.uniform(0.90,1.1)}
                                    }
        modulation_ACh = {'intr':   {'naf': np.random.uniform(1.0,1.2),
                                    'kir': np.random.uniform(0.5,0.7),
                                    'can': np.random.uniform(0.65,0.85),
                                    'cal12': np.random.uniform(0.3,0.7),
                                    'cal13': np.random.uniform(0.3,0.7),
                                    'Im': np.random.uniform(0.0,0.4)
                                    },
                         'syn':    {'NMDA':np.random.uniform(1.0,1.05), 
                                    'AMPA':np.random.uniform(0.99,1.01), 
                                    'GABA':np.random.uniform(0.99,1.01)} 
                                    }
        
        res[i][ci] = {'factors':{'da': modulation_DA, 'ach': modulation_ACh }  }
    
        # ACh factor sets used == 0 and 2 (since giving less inward rectification and higher excitability)
        V = {}
        for key in modulation_ACh['intr']:
            V[key] = transient
        for key in ['glut','gaba']:
            V[key] = transient
        
        # Master seed -> reset for all factor combinations
        #np.random.seed(seed=1000) # TODO no seed... will bg be different over workers?
        
        # set bg (constant+ramping). 
        fgaba=24.0
        Syn, nc, ns, mapper = use.set_ramping_stimuli(cell, [],low=modOnTime,high=modOnTime+tau)
        RAMP                = {'s':Syn.copy(), 'nc':nc.copy(), 'ns':ns.copy}
        Syn, nc, ns         = use.set_bg_noise( cell, fglut=12.0, fgaba=fgaba)
        
        # Clean
        DA = modulate.DA(   cell, modulation_DA['intr'],
                            modulation='uniform',
                            play=V,
                            syn_dict = modulation_DA['syn']
                            )
        ACh = modulate.ACh( cell, modulation_ACh['intr'],
                            shift_kaf=0,
                            play=V,
                            syn_dict = modulation_ACh['syn']
                            )
        
        for bg in range(nbg):    
            print(cell_index, c, bg)
            # finalize and run
            h.finitialize(-80)
            while h.t < tstop:
                h.fadvance()
            # downsample data and store in array
            if not 'time' in res:
                res['time'] = [t for ind,t in enumerate(tm) if ind%4 == 0]
            res[i][bg] = [vm[ind]   for ind,t in enumerate(tm) if ind%4 == 0]
        
        
    # save
    with open('inVivo_ramping_{}_{}_model{}.json'.format(conditions[ci],tag,cell_index), 'wt') as f:
        json.dump(res,f,indent=4)


# if run from terminal...   ===============================================================
if __name__ == "__main__":
    
    # select which model to use
    cell_index  =   rank%34
    cond_id     =   0 # ACh+DA
    
    run_model(cell_index, cond_id)
    
    


    
