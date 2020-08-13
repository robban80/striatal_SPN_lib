
# Minimal example on how to use the model libraries used in:
#   Lindroos & Hellgren Kotaleski 2020


from   neuron           import h
import MSN_builder          as build
import pickle
import matplotlib.pyplot    as plt

# Load model mechanisms
import neuron               as nrn
nrn.load_mechanisms('mechanisms/single')

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

# specs
specs = {'dspn': {
                    'N': 71,
                    'lib': 'Libraries/D1_71bestFit_updRheob.pkl',
                    'par': 'params_dMSN.json',
                    'morph': 'Morphologies/WT-dMSN_P270-20_1.02_SGA1-m24.swc'},
         'ispn': {
                    'N': 34,
                    'lib': 'Libraries/D2_34bestFit_updRheob.pkl',
                    'par': 'params_iMSN.json',
                    'morph': 'Morphologies/WT-iMSN_P270-09_1.01_SGA2-m1.swc'}
        }
        




# chose cell type ('ispn' or 'dspn') and model id(s) to simulate
#---------------------------------------------------------------
cell_type         = 'dspn'    # 'dspn'/'ispn'
model_iterator    = range(5)  # range(specs[cell_type]['N']) gives all models





# open library (channel distributions etc)   
with open(specs[cell_type]['lib'], 'rb') as f:
    model_sets = pickle.load(f, encoding="latin1")

# simulate model(s)
OUT = {}
for cell_index in model_iterator: 
    
    # initiate cell
    cell = build.MSN(  params=specs[cell_type]['par'],
                       morphology=specs[cell_type]['morph'],
                       variables=model_sets[cell_index]['variables']   )
    
    
    # set current injection
    rheobase        =   model_sets[cell_index]['rheobase']
    Istim           =   h.IClamp(0.5, sec=cell.soma)
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
    

# plot
for cell_index in OUT:       
    plt.plot(OUT[cell_index]['tm'], OUT[cell_index]['vm'], \
        label='mdl:{} rhb:{:.0f}'.format(cell_index,OUT[cell_index]['rheo']))
plt.legend()
plt.show()
        

