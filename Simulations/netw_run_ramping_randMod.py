
# script used for simulating in vitro activation of clustered synaptic input
# This script is based on plateaus_in_vivo_saveFeatures.py (in complex spikes repo)

from neuron import h

import sys
sys.path.insert(0, '../')
        
import json, pickle, glob
import numpy                as np
import CELL_builder_netw    as build
import common_functions     as use
import modulation_lib       as modulate
from   datetime             import date

# Load model mechanisms
#import neuron               as nrn
#nrn.load_mechanisms('../mechanisms/network/')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')


cells_dirs =    {
                'dspn': [
                    'str-dspn-e150917_c10_D1-mWT-P270-20-v20190521',
                    'str-dspn-e150917_c6_D1-m21-6-DE-v20190503',
                    'str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508',
                    'str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521'
                    ],
                'ispn': [
                    'str-ispn-e150917_c11_D2-mWT-MSN1-v20190603',
                    'str-ispn-e160118_c10_D2-m46-3-DE-v20190529',
                    'str-ispn-e150908_c4_D2-m51-5-DE-v20190611',
                    'str-ispn-e151123_c1_D2-mWT-P270-09-v20190527'
                    ]
                }


pc = h.ParallelContext()
id = int(pc.id())

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# import matplotlib.pyplot as plt

def write2file(text):
    with open("Outfile.txt", "a") as text_file:
        print(text, file=text_file)


def draw_random_factor(ct):
    if ct == 'ispn':
        modulation_DA = {'intr':   {'naf_ms': np.random.uniform(0.95,1.1),
                                    'kaf_ms': np.random.uniform(1.0,1.1),
                                    'kas_ms': np.random.uniform(1.0,1.1),
                                    'kir_ms': np.random.uniform(0.8,1.0),
                                    'can_ms': np.random.uniform(0.9,1.0),
                                    'car_ms': np.random.uniform(0.6,0.8),
                                    'cal12_ms': np.random.uniform(0.7,0.8),
                                    'cal13_ms': np.random.uniform(0.7,0.8)
                                    },
                         'syn':    {'NMDA':np.random.uniform(0.85,1.05), 
                                    'AMPA':np.random.uniform(0.7,0.9), 
                                    'GABA':np.random.uniform(0.90,1.1)}
                                    }
        modulation_ACh = {'intr':   {'naf_ms': np.random.uniform(1.0,1.2),
                                    'kir_ms': np.random.uniform(0.5,0.7),
                                    'can_ms': np.random.uniform(0.65,0.85),
                                    'cal12_ms': np.random.uniform(0.3,0.7),
                                    'cal13_ms': np.random.uniform(0.3,0.7)
                                    },
                         'syn':    {'NMDA':np.random.uniform(1.0,1.05), 
                                    'AMPA':np.random.uniform(0.99,1.01), 
                                    'GABA':np.random.uniform(0.99,1.01)} 
                                    }    
    elif ct == 'dspn':
        modulation_DA = {'intr':   {'naf_ms': np.random.uniform(0.6,0.8),
                                    'kaf_ms': np.random.uniform(0.75,0.85),
                                    'kas_ms': np.random.uniform(0.65,0.85),
                                    'kir_ms': np.random.uniform(0.85,1.25),
                                    'can_ms': np.random.uniform(0.2,1.0),
                                    'cal12_ms': np.random.uniform(1.0,2.0),
                                    'cal13_ms': np.random.uniform(1.0,2.0)
                                    },
                         'syn':    {'NMDA':np.random.uniform(1.2,1.4), 
                                    'AMPA':np.random.uniform(1.0,1.3), 
                                    'GABA':np.random.uniform(0.80,1.2)}
                                    }
        modulation_ACh = {'intr':   {'naf_ms': np.random.uniform(1.0,1.2),
                                    'kir_ms': np.random.uniform(0.8,1.0),
                                    'can_ms': np.random.uniform(0.65,0.85),
                                    'cal12_ms': np.random.uniform(0.3,0.7),
                                    'cal13_ms': np.random.uniform(0.3,0.7)
                                    },
                         'syn':    {'NMDA':np.random.uniform(1.0,1.05), 
                                    'AMPA':np.random.uniform(0.99,1.01), 
                                    'GABA':np.random.uniform(0.99,1.01)} 
                                    }    
    else:
        raise('modulation factors for type: {}'.format(ct))
    
    return modulation_DA, modulation_ACh

def run_model(  cell_type=None,
                mdl_ID=0,
                ci=0         ): 
    
    modeldir    =   '../Striatal_network_models/{}/{}'.format(cell_type, cells_dirs[cell_type][mdl_ID])
    
    par         =   '{}/parameters_with_modulation.json'.format(modeldir)
    mech        =   '{}/mechanisms.json'.format(modeldir)
    protocols   =   '{}/protocols.json'.format(modeldir)
    morphology  =   glob.glob(modeldir+'/*.swc')[0] # ONLY 1 swc file / model allowed.
    
        
    # initiate cell
    cell = build.CELL(  params=par,
                        mechanisms=mech,
                        morphology=morphology,
                        replace_axon=True  )    

    
    # 1s low input; 0.3s ramping input; 0.2s high input
    tstop       = 1500
    modOnTime   = 1000
    nbg         = 5
    tau         = 300
    conditions  = ['ACh+DA','ACh','DA','ctrl']
    c           = conditions[ci]
    
    # set current injection
    with open(protocols) as file:
        prot = json.load(file)
    
    # select first spiking prot
    key    = sorted(prot.keys())[0]
    st     = prot[key]['stimuli'][1]
    stim   =   h.IClamp(0.5, sec=cell.soma)
    stim.amp    =   st['amp']
    stim.delay  =   0
    stim.dur    =   tstop
    
    
    # record vectors
    tm  = h.Vector()
    vm  = h.Vector()
    tm.record(h._ref_t)
    vm.record(cell.soma(0.5)._ref_v)
    
    # transient to play
    transient = h.Vector([use.alpha(ht, modOnTime, 1, tau) if ht >= modOnTime else 0 for ht in np.arange(0,tstop,h.dt)])
    
    
    # result dict
    res = {}
    today = date.today() 
    tag = np.random.randint(9999)
    
    for i in range(20):
    
        # draw random factors # TODO add Im to models and ACh modulation
        modulation_DA, modulation_ACh = draw_random_factor(cell_type)
        
        res[i] = {'par': {'factors':{'da': modulation_DA, 'ach': modulation_ACh}, 'date':today.strftime("%b-%d-%Y"), 'cond':c }  }
    
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
        if 'DA' in c:
            DA = modulate.DA(   cell, modulation_DA['intr'],
                                modulation='uniform',
                                play=V,
                                syn_dict = modulation_DA['syn']
                                )
        if 'ACh' in c:
            ACh = modulate.ACh( cell, modulation_ACh['intr'],
                                shift_kaf=0,
                                play=V,
                                syn_dict = modulation_ACh['syn']
                                )
                            
        
        for bg in range(nbg):    
            print(cell_type, c, i, bg)
            # finalize and run
            h.finitialize(-80)
            while h.t < tstop:
                h.fadvance()
            # downsample data and store in array
            if not 'time' in res:
                res['time'] = [t for ind,t in enumerate(tm) if ind%4 == 0]
            res[i][bg] = [vm[ind]   for ind,t in enumerate(tm) if ind%4 == 0]
        
        # save
        with open('netw_{}_ramping_{}_{}_model{}.json'.format(cell_type, conditions[ci],tag,mdl_ID), 'wt') as f:
            json.dump(res,f,indent=4)
     
         
        
        
# if run from terminal...   ===============================================================
if __name__ == "__main__":
    
    # select which model to use
    cell_type   = 'dspn' # ['dspn','ispn']
    mdl_ID      = rank%4
    # modulation paradigm (0-3; ACh+DA, ACh, DA, ctrl)
    ci          =   0 # ? int(np.floor(id/40))
    
    run_model(cell_type, mdl_ID, ci)
    
    
    
    
