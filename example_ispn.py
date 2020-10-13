
# example script for showing how to use the model library

'''
select models to simulate below (mdls). 
    By default the five first will be simulated with a current injection
    of rheobase + random(0-60) pA
'''


from   neuron           import h
import MSN_builder          as build
import pickle
import matplotlib.pyplot    as plt
import numpy                as np

# Load model mechanisms
import neuron               as nrn
nrn.load_mechanisms('mechanisms/single')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')



# open channel distribution parameter library    
with open('Libraries/D2_34bestFit_updRheob.pkl', 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1") 


# basic parameters and morphology
par         =   './params_iMSN.json'
morphology  =   './Morphologies/WT-iMSN_P270-09_1.01_SGA2-m1.swc'

# select model index here (use range or a list, e.g. [0,1,4,20] or [0])
mdls        =   range(5)

OUT = {}
for cell_index in mdls:
    
    parameters      =   model_sets[cell_index]['variables'] 
    rheobase        =   model_sets[cell_index]['rheobase']

    # initiate cell -
    cell = build.MSN(  params=par,                  \
                       morphology=morphology,       \
                       variables=parameters         )
    
    
    # set current injection
    Istim           = h.IClamp(0.5, sec=cell.soma)
    Istim.delay     =   100
    Istim.dur       =   1000
    randint         =   np.random.randint(60)
    Istim.amp       =   (rheobase+randint) *1e-3
    
    # record vectors
    tm  = h.Vector()
    tm.record(h._ref_t)
    vm  = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)
    
    # run simulation
    h.finitialize(-80)
    
    # run simulation
    while h.t < 1000:
        h.fadvance()
    
    OUT[cell_index] = {'tm':tm.to_python(), 'vm':vm.to_python, 'rheo':rheobase}
    

# plot
for cell_index in OUT:       
    plt.plot(OUT[cell_index]['tm'], OUT[cell_index]['vm'], \
        label='mdl:{} rhb:{:.0f}'.format(cell_index,OUT[cell_index]['rheo']))
plt.legend()
plt.show()
        

