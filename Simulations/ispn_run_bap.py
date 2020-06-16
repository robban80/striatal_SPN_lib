
# run FI curve stim for all ispn's. rheobase -> rheobase+90 (dI = 30 pA)


from   neuron           import h

import sys
sys.path.insert(0, '../')

import MSN_builder          as build
import common_functions     as use
import numpy                as np
import pickle, json


# Load model mechanisms
import neuron               as nrn
nrn.load_mechanisms('../mechanisms/single')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')


# curve_fit function
from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a * np.exp( (x-b) / c )



# open channel distribution parameter library    
with open('../Libraries/D2_34bestFit_updRheob.pkl', 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1") 


# basic parameters and morphology
par         =   '../params_iMSN.json'
morphology  =   '../Morphologies/WT-iMSN_P270-09_1.01_SGA2-m1.swc'

N       = 34
res     = {}

for cell_index in range(N):
    
    parameters      =   model_sets[cell_index]['variables'] 
    rheobase        =   model_sets[cell_index]['rheobase']
    
    # initiate cell -
    cell = build.MSN(  params=par,                  \
                       morphology=morphology,       \
                       variables=parameters         )
    
    ca_rec = {}
    for sec in cell.dendlist:
        for seg in sec:
            d = h.distance(seg.x, sec=sec)
            ca_rec[seg] = {'cal':h.Vector(),'ca':h.Vector(),'dist':d}
            ca_rec[seg]['cal'].record(seg._ref_cali)
            ca_rec[seg]['ca'].record(seg._ref_cai)
    # TODO test above record implementation
    # - implement ca averaging and save to file
    
    # set current injection
    Istim           = h.IClamp(0.5, sec=cell.soma)
    Istim.delay     =   100
    Istim.dur       =   2
    Istim.amp       =   2 # 2000 pA = 2 nA
    
    # record vectors
    tm  = h.Vector()
    tm.record(h._ref_t)
    vm  = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)
    
    # run simulation
    h.finitialize(-80)
    
    # run simulation
    while h.t < 200:
        h.fadvance()
    
    # check if spike is triggered:
    if max(vm) < 20:
        import matplotlib.pyplot as plt
        plt.plot(tm,vm)
        plt.plot(tm[3999],vm[3999], 'or', ms=8)
        plt.show()
        gggh
    
    # make into lists
    caList   = np.zeros(len(ca_rec))
    for k,key in enumerate(ca_rec): # (max(cai)+max(cali)) - (cai[0]+cali[0])
        caList[k]   = (max(ca_rec[key]['ca'].to_python()[3999:])+max(ca_rec[key]['cal'].to_python()[3999:])) \
                    - (ca_rec[key]['ca'][3999]+ca_rec[key]['cal'][3999])
    
    res[cell_index] = list(caList)
    if 'dist' not in res:
        distList = np.zeros(len(ca_rec))
        for k,key in enumerate(ca_rec):
            distList[k] = ca_rec[key]['dist']
        res['dist'] = list(distList)
        '''
        # normalize to random (first) compartment within 30-40 um somatic dist
        maxdist     = 210
        index       = np.argsort( distList )
        norm_ind,D  = next(i for i in enumerate(distList) if i[1] < 40 and i[1] > 30) 
        dist        = [distList[i] for i in index \
                                        if distList[i] >= D and distList[i] < maxdist]
        res['dist'] = list(dist)'''
    
    '''   
    Csum        = [caList[i]/caList[norm_ind] for i in index \
                                        if distList[i] >= D and distList[i] < maxdist]
    # get regression lines
    popt, pcov = curve_fit( func, 
                            dist, 
                            Csum, 
                            p0=[1,40,-15])
    
    # save
    res[model_version] = list(func(dist, *popt))
    '''
    
res['tm'] = tm.to_python()
res['vm'] = vm.to_python()

outFile = '../Validation_data/ispn_extracted_bap.json'
#serialize_and_save_json(res, outFile)
with open(outFile, 'wt') as f:
    json.dump(res,f,indent=4)
    

        

