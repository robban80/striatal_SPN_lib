#
'''
Example script showing how to run network models outside of BluePyOpt. 
    Models are optained from Hjorth et al., 2020 (PNAS)
    
--------------------------------------------------------------------------------------

The models use uniform channel distribution in the dendrites due to earlier bug.


Implementation of conversion done by:
    Robert Lindroos (RL) <robert.lindroos at ki.se>
 
Optimized models and original framework by:  
    Alexander Kozlov (AK) <akozlov at kth.se>

'''

from __future__ import print_function, division


from neuron import h
import glob, json, pickle
import numpy                as np
import CELL_builder_netw    as build
#import common_functions     as use

# Load model mechanisms
import neuron               as nrn
nrn.load_mechanisms('./mechanisms/network/')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')


# TODO: add chin and lts?
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
    
        
# main ==================================================================================

def main(       cell_type=None,
                mdl_ID=0         ): 
    
    modeldir    =   './Striatal_network_models/{}/{}'.format(cell_type, cells_dirs[cell_type][mdl_ID])
    
    par         =   '{}/parameters_with_modulation.json'.format(modeldir)
    mech        =   '{}/mechanisms.json'.format(modeldir)
    protocols   =   '{}/protocols.json'.format(modeldir)
    morphology  =   glob.glob(modeldir+'/*.swc')[0] # ONLY 1 swc file / model allowed.
    
        
    # initiate cell
    cell = build.CELL(  params=par,
                        mechanisms=mech,
                        morphology=morphology,
                        replace_axon=True  )
                
    # THIS PART IS OPTIONAL oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # set input here
    
    # set current injection
    with open(protocols) as file:
        prot = json.load(file)
    
    # select first spiking prot
    all_keys    = sorted(prot.keys())
    key         = all_keys[0]
    
    i=1
    while 'sub' in key:
        key = all_keys[i]
        i += 1
    print(key)
    
    stim        =   h.IClamp(0.5, sec=cell.soma)
    s0          =   h.IClamp(0.5, sec=cell.soma)
    
    for stim_prot, stimuli, j in zip(prot[key]['stimuli'], [stim,s0], [0,1]): 
        stimuli.amp    =   stim_prot['amp']
        stimuli.delay  =   [200,0][j]
        stimuli.dur    =   stim_prot['duration']
    
    
    # oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
     
    # record vectors: set recordings here
    tm  = h.Vector()
    tm.record(h._ref_t)
    vm  = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)
    
    tstop  = 600    # sim time (ms)
    
    h.finitialize(cell.v_init)
    
    # run simulation
    while h.t < tstop:
        h.fadvance()
        
    # save output ------------------------------------------------------------------------
    time    = tm.to_python()
    voltage = vm.to_python()
    
    return [time, voltage]
    
        

# Start the simulation and save results
if __name__ == "__main__":
    
    res = {}
    
    cell_types  = ['dspn','ispn']
    
    # model id's. 4 for each spn
    for cell_type in cell_types:
        res[cell_type] = {}
        for mdl_ID in range(4):
        
             res[cell_type][mdl_ID] = main( mdl_ID=mdl_ID,
                                            cell_type=cell_type )  
    
    # save traces
    save_dict = {'time': res[cell_type][mdl_ID][0]}
    for key,item in res.items():
        save_dict[key] = {}
        for k,trace in item.items():
            save_dict[key][k] = trace[1]
    
    # UNCOMMENT TO SAVE
    #with open('val_data.json', 'w') as outfile:
    #    json.dump(save_dict, outfile)
    
    # plot
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(2,1)
    for i,cell_type in enumerate(cell_types):
        for mdl_ID in range(4):        
            t = res[cell_type][mdl_ID][0]
            v = res[cell_type][mdl_ID][1]
            ax[i].plot(t,v, 'k')
        ax[i].set_ylim([-90,40])
    plt.show() 
    
        
        
    
                                                    
    
                                                    
                                                    
                                                    
                                                    
    
    
    
          
    
        

