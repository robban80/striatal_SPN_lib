
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

# import matplotlib.pyplot as plt

def write2file(text):
    with open("Outfile.txt", "a") as text_file:
        print(text, file=text_file)


def run_model(  cond,
                pattern,
                cell_type=None,
                mdl_ID=0         ): 
    
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
    modfunctype = 'sig'
    
    mod_nmda    = 0
    sf_nmda     = 1.0
    
    da=0; ach=0
    
    print(cell_type, mdl_ID)
    
    # channel distribution parameter library   
    '''
    path            = 'Cholinergic_modulation'
    
    # open factors
    with open('Dopamine_modulation/dspn_DA_factors_v2.json', 'r') as f:
        DAfactors   = json.load(f)
    with open('{}/dspn_ACh_factors_v2.json'.format(path), 'r') as f:
        ACHfactors  = json.load(f)
    '''
    # load random sets
    mixActPat   = use.load_obj('../Libraries/mixActPatRestricted16.pkl')
    # ====================================================================================
                    
    # set stimuli (including modulation)
    rand, ncon, stim, spike_list = use.set_mixed_stimuli(   cell,               \
                                                            mixActPat,          \
                                                            pattern,            \
                                                            ISI=1               )
    
    # record vectors
    tm  = h.Vector()
    vm  = h.Vector()
    tm.record(h._ref_t)
    vm.record(cell.soma(0.5)._ref_v)
    
    # transient to play
    '''
    if modfunctype == 'sig':
        transient = h.Vector([use.sigmoid(ht, modOnTime, slope=-5) for ht in np.arange(0,tstop,h.dt)])
    elif modfunctype == 'alpha':
        transient = h.Vector([use.alpha(ht, modOnTime, slope=-5) for ht in np.arange(0,tstop,h.dt)])
    elif modfunctype == 'step':
        transient = h.Vector([1 if ht > modOnTime else 0 for ht in np.arange(0,tstop,h.dt)])
    else: raise Exception('Error, no modulation function "{}"'.format(modfunctype))
    '''
    
    VAL = {}
       
    for i in range( 1+len(mixActPat['steps']) ):
        
        # (re-)set bg and shift annealing
        if i > 0:
            try:
                spike_list, stim, ncon  = use.mixed_stimuli_annealing(mixActPat, pattern, i-1, stim, ncon, rand, spike_list)
            except:
                break
            if i%2 == 1: continue   
        
        VAL[i] = {}
        
        
        
        # set modulation =================================================================
        # modulation factors
        MOD = False
        if MOD:
            modulation_ACh = {}; V = {}
            shift_kaf = 0; kaff = 20    # base values. automatically updated if in factors
            for key,val in ACHfactors[str(ach)].items():
                if key == 'kaf': shift_kaf = 1; kaff = ACHfactors[str(i)]['kaf']
                else: modulation_ACh[key] = val
                V[key] = transient
            modulation_DA = {}
            for key,val in DAfactors[str(da)].items():
                modulation_DA[key] = val
                V[key] = transient
            for key in ['glut','gaba']:
                V[key] = transient   
            
            # Clean
            if 'ACh' in locals(): ACh._reset_mod()
            if 'DA'  in locals(): DA._reset_mod()
            # Set modulation
            if cond == 'DA':
                DA = modulate.DA(   cell, modulation_DA,
                                    play=V,
                                    syn_dict={'NMDA':1.3, 'AMPA':1.2, 'GABA':0.8}
                                    )
            elif cond == 'ACh':
                ACh = modulate.ACh( cell, modulation_ACh,
                                    shift_kaf=shift_kaf,
                                    mv_shift_kaf=kaff,
                                    play=V,
                                    syn_dict={'NMDA':sf_nmda, 'AMPA':1.0, 'GABA':1.0}
                                    )
            elif cond == 'ACh+DA':
                DA = modulate.DA(   cell, modulation_DA,
                                    play=V,
                                    syn_dict={'NMDA':1.3, 'AMPA':1.2, 'GABA':0.8}
                                    )
                ACh = modulate.ACh( cell, modulation_ACh,
                                    shift_kaf=shift_kaf,
                                    mv_shift_kaf=kaff,
                                    play=V,
                                    syn_dict={'NMDA':sf_nmda, 'AMPA':1.0, 'GABA':1.0}
                                    )
            elif cond == 'ctrl': pass
            else: raise Exception('Error, no condition "{}"'.format(cond))
            
        # ================================================================================
        
        # finalize and run ----------------------------------------------------------------
        h.finitialize(cell.v_init)
        while h.t < tstop:
            h.fadvance()   
            
        skiptT = 900     
        # shift spike traces and skip first 500 ms
        t   = [x-1000    for x in tm            if x    >= skiptT]
        v   = [vm[x[0]]  for x in enumerate(tm) if x[1] >= skiptT]
        
        VAL[i] = { 'vm': v }
    VAL['time'] = t
    
    return VAL
     
         
        
        
# if run from terminal...   ===============================================================
if __name__ == "__main__":
    
    N = id%40 #40
    
    # select which model to use
    ci          =   int(np.floor(id/40)) #int(h.modind) #int(np.floor(id/54))
    conditions  =   ['ctrl','ACh','DA','ACh+DA']
    cond        =   conditions[ci]
    
    today = date.today() 
    RES = {'par':    {  'date':today.strftime("%b-%d-%Y")}}
    
    cell_types  = ['dspn','ispn']
    
    for cell_type in cell_types:
        RES[cell_type] = {}
        for mdl_ID in range(4):                 # model id's. 4 for each spn
            RES[cell_type][mdl_ID] = {}
            for pattern in range(N,N+1):
                VAL = run_model(cond, pattern, cell_type, mdl_ID)
                #write2file('model: {}, pattern: {}'.format(cell_index, pattern))
                RES[cell_type][mdl_ID][pattern] = VAL
    
    # save traces
    save_dict = {'time': RES['dspn'][0][0][0]}
    for key,item in RES.items():
        if key == 'par':continue
        save_dict[key] = {}
        for k2,i2 in item.items():
            save_dict[key][k2] = {}
            print(k2)
            for k3,i3 in i2.items():
                print('-',k3)
                save_dict[key][k2][k3] = i3
    
    # UNCOMMENT TO SAVE
    with open('inVitro_mixedPattern_res.json', 'w') as outfile:
        json.dump(save_dict, outfile)
    '''            
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(2,4)
    for i,cell_type in enumerate(cell_types):
        for mdl_ID in range(4):
            for pattern in range(N):
                VAL = RES[cell_type][mdl_ID][pattern]
                for j in VAL:
                    if j == 'time': continue
                    ax[i,mdl_ID].plot(VAL['time'], VAL[j]['vm'])
    plt.show()'''
        
    
    #use.save_obj(RES, 'inVitro2020_mixedPattern16_D1_V{}_mod{}'.format(cell_index,cond) )
    
    
