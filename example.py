
# script illustrating the use the model library (dspn)


from   neuron           import h
import MSN_builder          as build
import common_functions     as fun
import pickle
import matplotlib.pyplot    as plt

# Load model mechanisms
import neuron               as nrn
nrn.load_mechanisms('mechanisms_with_pointer')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')



# open channel distribution parameter library    
with open('Libraries/D1_71bestFit.pkl', 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1") 


# basic parameters and morphology
par         =   './params_dMSN.json'
morphology  =   'morphologies/WT-dMSN_P270-20_1.02_SGA1-m24.swc'


# chose model(s) here. 
#-----------------------------------------
OUT = {}
for cell_index in [2]:#range(71):
    
    parameters      =   model_sets[cell_index]['variables'] 
    rheobase        =   model_sets[cell_index]['rheobase']

    # initiate cell -
    cell = build.MSN(  params=par,                  \
                       morphology=morphology,       \
                       variables=parameters         )
    
    h.topology()
    mod_list = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can' ]
    
    casc    =   h.D1_reduced_cascade2_0(0.5, sec=cell.soma)
    pointer =   casc._ref_Target1p
    fun.set_pointers(cell, pointer, mod_list)
   
    # set current injection
    Istim           = h.IClamp(0.5, sec=cell.soma)
    Istim.delay     =   100
    Istim.dur       =   1000
    Istim.amp       =   (rheobase) *1e-3
    
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
    

for cell_index in OUT:       
    plt.plot(OUT[cell_index]['tm'],OUT[cell_index]['vm'], label='%-%'.format(cell_index,rheobase))
    #plt.title('model v: %d; rheobase: %d' % (cell_index, rheobase))
    #plt.savefig('model_%02d.png' %(cell_index), format='png')
plt.legend()
plt.show()
        

