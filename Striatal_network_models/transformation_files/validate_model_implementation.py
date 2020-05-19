#
'''
Transformation and random retuning of models from striatal network.
    Based on older script used to randomly alter the model from Lindroos et al., 2018
    
------------------------------------------------------------------------------------------
'''
from __future__ import print_function, division
# TODO: udpate below
'''
Master script used for random tuning of dMSN and iMSN. 

-The tuning is done by randomly chosing dendritic distribution parameters for the 
following channels:

    ['sk', 'can' 'cat32', 'cat33', 'kir', 'kas', 'kaf', 'naf']

from a sigmoidal distribution (except for kir and sk that are uniform):

    y = a4 + a5/(1 + np.exp((x-a6)/a7) )

kaf:
a4  1
a5  np.random.uniform(0.1,0.9)
a6  np.random.uniform(1,130)
a7  np.random.uniform(-3,-70)

naf:
a4  1-a5
a5  np.random.uniform(0.0,1.0)
a6  np.random.uniform(20,50)
a7  np.random.uniform(1,50)

and then validating the somatic excitability to Planert et al 2013 (FI curve)
and dendritic excitability to Day et al 2008 (bAP induced Ca change).

-Cell type is randomly chosen. 
The differences between the cell type protocoles are:

* Kaf channel "base value" (0.11 vs 0.06 in dMSN and iMSN, respectively)
* morphology   

Robert Lindroos (RL) <robert.lindroos at ki.se>
 
Based on original Neuron implementation by  
Alexander Kozlov <akozlov at kth.se>
'''



from neuron import h
import glob, json, pickle
import numpy                as np
import CELL_builder         as build
#import common_functions     as use

# Load model mechanisms
import neuron               as nrn
nrn.load_mechanisms('../../mechanisms/network/') # TODO update this path

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')


# global result dicts
RES             = {}
VARIABLES       = {}
CA              = {}

# TODO
'''
-access models from respective model catalog
-add channel distribution from param file
-test
-brute force alter these models and add to lib 
-(make 1 lib?)
'''

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
    
    modeldir    =   '../{}/{}'.format(cell_type, cells_dirs[cell_type][mdl_ID])
    
    par         =   '{}/parameters_with_modulation.json'.format(modeldir)
    mech        =   '{}/mechanisms.json'.format(modeldir)
    protocols   =   '{}/protocols.json'.format(modeldir)
    morphology  =   glob.glob(modeldir+'/*.swc')[0]
    
        
    # initiate cell
    cell = build.CELL(  params=par,
                        mechanisms=mech,
                        morphology=morphology,
                        replace_axon=True  )
                
    print(modeldir.split('/')[-1])
    #for sec in cell.axonlist:
    #    h.psection(sec=sec)
    #h.psection(sec=cell.soma)
    #print(cell.soma(0.5).area())
    '''
    nsec = 0
    nseg = 0
    for sec in cell.axonlist:
        nsec += 1
        nseg += sec.nseg
    print('AXON: sections: {}, segments {}'.format(nsec, nseg))
    nsec = 1
    nseg = cell.soma.nseg
    print('SOMA: sections: {}, segments {}'.format(nsec, nseg))
    for sec in cell.dendlist:
        nsec += 1
        nseg += sec.nseg
    print('DEND: sections: {}, segments {}'.format(nsec, nseg))
    for i,sec in enumerate(cell.dendlist):
        for j,seg in enumerate(sec):
            for mech in seg:
                if mech.name() not in ['naf_ms']: continue
                try:
                    print(i,j, mech.name(), mech.gbar, seg.area())
                except:
                    pass    '''
                       
    
    
    
    # set current injection
    with open(protocols) as file:
        prot = json.load(file)
    
    # select first spiking prot
    all_keys    = sorted(prot.keys())
    i           = 0
    key         = all_keys[i]
    
    while 'sub' in key:
        i  += 1
        key = all_keys[i]
    print(key)
    stim        =   h.IClamp(0.5, sec=cell.soma)
    s0          =   h.IClamp(0.5, sec=cell.soma)
    
    for stim_prot, stimuli, j in zip(prot[key]['stimuli'], [stim,s0], [0,1]): 
        stimuli.amp    =   stim_prot['amp']
        stimuli.delay  =   [200,0][j]
        stimuli.dur    =   stim_prot['duration']
    
    # record vectors
    tm  = h.Vector()
    tm.record(h._ref_t)
    vm  = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)
    
    tstop       = 600
    
    
              
    # solver------------------------------------------------------------------------------            
    cvode = h.CVode()
    
    h.finitialize(cell.v_init)
    
    # run simulation
    while h.t < tstop:
                
        h.fadvance()
        
    # save output ------------------------------------------------------------------------
    
    time = tm.to_python()
    voltage = vm.to_python()
    
    return [time, voltage]
    
        

# Start the simulation and save results
# Function needed for HBP compability  ===================================================
if __name__ == "__main__":
    
    res = {}
    
    cell_types  = ['dspn','ispn']
    
    # model id's. 4 for each spn
    for cell_type in cell_types:
        res[cell_type] = {}
        for mdl_ID in range(4):
        
             res[cell_type][mdl_ID] = main( mdl_ID=mdl_ID,
                                            cell_type=cell_type )  
    
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots(2,1)
    
    path2control = '../../../Alex_model_repo/models/optim/Dopamine/Analysis/Results/'
    with open('{}dspn_res_org.pkl'.format(path2control), 'rb') as f:
        controlD = pickle.load(f)
    with open('{}ispn_res_org.pkl'.format(path2control), 'rb') as f:
        controlI = pickle.load(f)

    ctrl = {'ispn':controlI, 'dspn':controlD}

    mapid = {'dspn':[3,1,0,2] , 'ispn':[1,3,0,2] }
    
    # save traces
    save_dict = {'time': res[cell_type][mdl_ID][0]}
    for key,item in res.items():
        save_dict[key] = {}
        for k,trace in item.items():
            save_dict[key][k] = trace[1]
    
    with open('val_data.json', 'w') as outfile:
        json.dump(save_dict, outfile)
    
    for i,cell_type in enumerate(cell_types):
        for mdl_ID in range(4):        
            t = res[cell_type][mdl_ID][0]
            v = res[cell_type][mdl_ID][1]
            ax[i].plot(t,v, 'k')
            data = list(ctrl[cell_type.lower()][0]['data'][mdl_ID]['control'].values())[0]
            ax[i].plot(data['t'],data['v'] , c='r', ls='--')
        ax[i].set_xlim([0,200])
        ax[i].set_ylim([-90,-80])
    plt.show() 
    
        
        
    
                                                    
    
                                                    
                                                    
                                                    
                                                    
    
    
    
          
    
        

