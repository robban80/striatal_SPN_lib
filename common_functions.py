'''
common functions for MSN model
'''


from    neuron  import          h
import  numpy                as np

import pickle
import json, codecs






def alpha(ht, tstart, gmax=1, tau=500):
    ''' 
    calc and returns a "magnitude" using an alpha function -> used for modulation 
        transients
    
    ht      = simulation time (h.t)
    tstart  = time when triggering the function
    gmax    = maximal amplitude of curve (default 1; transient must lie between 0-1)
    tau     = time constant of alpha function
    '''
    
    t   = (ht - tstart) / tau
    e   = np.exp(1-t)
    mag = gmax * t * e
    
    return mag


def sigmoid(ht, tstart, const=0, gmax=1, slope=-5):
    ''' 
    calc and returns a "magnitude" using an sigmoidal function -> used for modulation 
        transients. 
        mag = const + transient
    
    ht      = simulation time (h.t)
    tstart  = time when triggering the step
    const   = constant value (default 0; should be zero if used as modulation transient)
    gmax    = amplitude of transient (default 1; transient must lie between 0-1)
    slope   = slope of step (from constant to constant+transient).
                by setting a positive slope the function in stead will be
                    constant+transient -> constant   
    
    '''
    
    return const + gmax/(1+np.exp((ht-tstart)/slope))
   
    
def serialize_and_save_json(obj, name, mod=1):
    '''
    save dict including ndarray to json
    https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
    '''
    # also picke? Remove?
    with open(name.replace('json','pkl'), 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
    obj['ctrl'] = obj['ctrl'].tolist()
    if mod: obj['mod' ] = obj['mod' ].tolist()
    json.dump(  obj, 
                codecs.open(name, 'w', encoding='utf-8'), 
                separators=(',', ':'),  
                sort_keys=True, 
                indent=4)

def unserialize_from_json(name, mod=1):
    '''
    Reverse serialization of json into ndarray. To run on json files created using above
    https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
    ''',codecs
    obj_text = codecs.open(name, 'r', encoding='utf-8').read()
    loaded_data = json.loads(obj_text)
    loaded_data['ctrl'] = np.array(loaded_data['ctrl'])
    if mod:loaded_data['mod' ] = np.array(loaded_data['mod' ])
    return loaded_data
    

def save_vector(x, y, outfile):
    '''
    save vectors to file.
    
    x       = x-vector
    y       = y-vector
    outfile = file name to be used
    '''
    
    with open(outfile, "w") as out:
        for time, y in zip(x, y):
            out.write("%g %g\n" % (time, y))
            
            
            

def save_obj(obj, name ):
    '''
    functions used to save data in the pickle format. 
    Used to pickle dictionaries
    
    obj     = dictionary with data
    name    = name to be used, without file ending (.pkl will be added)
    '''
    
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        
        

def load_obj(name ):
    '''
    functions used to load data in the pickle format. 
    Used to un-pickle data from dictionaries
    
    name    = file name of the pickled data, including ending
    '''
    
    with open(name, 'rb') as f:
        return pickle.load(f) 
        
def load_json(name ):
    '''
    functions used to load data in the json format.
    
    name    = file name of the json data, including ending
    '''
    
    with open(name, 'rb') as f:
        return json.load(f)            
    
def check_fit(trace, cell='D1'):
    '''
    check_fit(trace, cell='D1')
    
    calculates the sum of squares error between model and experimental data from Day et al
    (2008), regaridng bAP induced calcium entry.
    
    trace = mean dCa as function of somatic distance
    cell  = msn subtype (D1; default, or D2) 
    '''
    
    if cell=='D1':
        # (dMSN)
        dist = np.arange(40,200,10)
        diff = np.exp( (dist-40)/-27.0 ) - trace
        
    else:
        # (iMSN)
        dist = np.arange(49,200,10)
        diff = np.exp( (dist-49)/-63.0 ) - trace
    
    fit = np.square(diff)
        
    return np.sum(fit)
    
    
        
    
def getSpikedata_x_y(x,y, threshold=0.0):
    ''' 
    There's probably a Neuron function for this (e.g. netcon)--use instead? 
    
    getSpikedata_x_y(x,y) -> return spike times as a list
    
    Extracts and returns the spikes from spike trace data.
    
    extraction algorithm:
    -extract list containing index for all points in y larger than threshold.
    -sorts out the index(es) that are the first one(s) crossing the threshold, i.e. the 
        first index of each spike, and stores the timing of the event (x[index]). 
    This is done by looping over all index and check if the index is equal to the previous 
        index + 1. 
    If not it is the first index of a spike.
        
    If no point is above threshold in the trace the function returns an empty list [].
    
    x           = time vector
    y           = vm vector
    threshold   = min amplitude used for spike detection.
                  If the amplitude of the spike is lower than threshold spikes won't be
                    detected.
                  If the Threshold is not crossed in the repolarizing phase, no new spike 
                    can be detected.
    '''
    
    #count = 0
    spikes = []
    
    # pick out index for all points above zero potential for potential trace
    spikeData = [i for i,v in enumerate(y) if v > threshold]

    # if no point above 0
    if len(spikeData) == 0:
        
        return spikes
	
    else:
        # pick first point of each individual transient (spike)...
        for j in range(0, len(spikeData)-1):
            if j==0:
                
                #count += 1
                spikes.append(x[spikeData[j]])
		
            # ...by checking above stated criteria
            elif not spikeData[j] == spikeData[j-1]+1:
                #count += 1
                spikes.append(x[spikeData[j]])
            
    return spikes 
    


def draw_random_variables(return_range=False):
    
    # def ranges
    
    # if return_ranges; return ranges
    if return_range:
        return 2
    
    else:
    
        # draw random variables
        
        # return variables
        return 1
        

        

def record_axial_current(node, exception=False):
    # creates record structure needed for calculation of axial current flowing between
    #         a node and its child compartments
    # node is a section
    # if exception in stem; skip
    
    AXIAL = {}      
    for stem in node.children(): 
        
        if exception:
            flag = False
            for string in exception:
                if string in stem.name():
                    flag = True
                    break
            if flag:
                continue
        
        AXIAL[stem.name()] = { 'stem':stem, 'ri':stem(1e-3).ri() }
        
        AXIAL[stem.name()]['vm'] = h.Vector()
        AXIAL[stem.name()]['vm'].record( stem(1e-3)._ref_v )
        
    return AXIAL
    

def plot_axial_current( AXIAL, tm, vm, ax, mark_section=False ):
    # Calculates and plots axial current from structure defined by  
    #       
    #       record_axial_current
    
    ALL = []
    for stem in AXIAL:
        
        diff = np.subtract( vm, AXIAL[stem]['vm'] )
        I    = np.divide( diff, AXIAL[stem]['ri'] )
        
        flag  = True
        color = 'b'
        if mark_section:
            
            # is this stem parent stem to section to mark?
            
            # create list of sections to plot
            tree = h.SectionList()
            tree.subtree(sec=AXIAL[stem]['stem'])
            
            for sec in tree:
                
                if sec.name() in mark_section:
                    ax.plot(tm, I, color='r', lw=1, alpha=0.3)
                    flag  = False
        
        if flag:
            ALL.append(I)
                
    s = np.zeros( len(I) )
    for a in ALL: 
        s = np.add( s, a )
        
    ax. plot(tm, s, color='b', alpha=0.3)
    
    
    
        
def FWHM(tm, Y, baseIndex):
    '''
    calculates full width half max of a voltage trace.
    The trace can only have one epsp, i.e. it is assumed that the curve is monotonic
    before and after the maximum.
    '''
    
    # get half-max vm
    half_max = (max(Y)+Y[baseIndex]) / 2.

    indx = []
    for i,y in enumerate(Y):
        
        # get all values above half max (assumes monotonic epsp in trace)
        if y - half_max >= 0:
            indx.append(i)
            
    return tm[indx[-1]] - tm[indx[0]]
    
    
    
def HMDur(tm, Y, baseIndex):
    '''
    calculates half max duration from end of stimulation of a voltage trace.
    -The trace can only have one epsp, i.e. it is assumed that the curve is monotonic
    before and after the maximum.
    -The end of the voltage stimulation is assumed to be after 140 ms (100 delay + 20x2 ms ISI)
    '''
    
    # get half-max vm
    half_max = (max(Y)+Y[baseIndex]) / 2.

    indx = []
    for i,y in enumerate(Y):
        
        # get all values above half max (assumes monotonic epsp in trace)
        if y - half_max >= 0:
            indx.append(i)
            
    return tm[indx[-1]] - 140



def Plateau_area(tm, Y, baseIndex):
    
    # set baseline to 0
    Y = np.subtract( Y, Y[baseIndex] ) 
    
    # return area under curve
    return np.trapz( Y[baseIndex:], dx=tm[1]-tm[0] )   
    



def random_synapse(ns, nc, Syn, sec, x,         \
                Type='glut',                    \
                NS_start=0,                     \
                NS_interval=1000.0/18.0,        \
                NS_noise=1.0,                   \
                NS_number=1000,                 \
                S_AN_ratio=1.0,                 \
                S_tau_dep=100,                  \
                S_U=1,                          \
                S_e=-60,                        \
                S_tau1=0.25,                    \
                S_tau2=3.75,                    \
                NC_delay=0,                     \
                NC_conductance=0.6e-3,          \
                NC_threshold=0.1,               \
                seed = None                     ):
    '''
    random_synapse(argument, *, **)
    
    ---arg n removed. used for setting multiple gaba synapses in same segment
    
    creates a synapse in the segment closest to x in section sec, and updates the dict
    containing the synapses (as well as netStim and NetCon objects).
    
    Use the Type argument to specify synapse mechanism:
        Type        mechanism       description
        glut        tmglut          glutamatergic (ampa+nmda) with short term depression (default) 
        gaba        gaba.mod        exp2syn with modulation function
        ampa        Exp2syn         deprived. No longer supported.
        gabaOld     Exp2syn         NEURON native synapse with beta type dynamics
        tmgabaa     tmgabaa         gabaergic with short term depression
    
    Any other Type than the above stated will result in an error.
        
    NS_arguments;       defines the NetStim object
    S_arguments;        defines the synapse mechanism
    NC_arguments;       defines the NetCon  object
    '''
    
    # create/set synapse in segment x of section
    if Type == 'tmglut':
        key                 = sec.name() + '_glut'
        Syn[key]            = h.tmGlut(x, sec=sec)
        Syn[key].nmda_ratio = S_AN_ratio
        Syn[key].tauR       = S_tau_dep
        Syn[key].U          = S_U
        
    elif Type == 'ampa':
        raise Exception('ampa in Exp2syn format no longer supported! In: common_functions -> random_synapse')
        
    elif Type == 'glut':
        key                 = sec.name() + '_glut'
        Syn[key]            = h.glutamate(x, sec=sec)
        Syn[key].ratio      = S_AN_ratio
        
    elif Type == 'gabaOld':
        key                 = sec.name() + '_gaba'  #+str(n)
        Syn[key]            = h.Exp2Syn(x, sec=sec)
        # what more?
        Syn[key].tau1       = S_tau1
        Syn[key].tau2       = S_tau2
        Syn[key].e          = S_e
    
    elif Type == 'gaba':
        # this will no mechanism error if gaba.mod not in among mechanisms.
        key                 = sec.name() + '_gaba'  #+str(n)
        Syn[key]            = h.gaba(x, sec=sec)
        
    elif Type == 'tmgabaa':
        key                 = sec.name() + '_gaba'
        Syn[key]            = h.tmGabaA(x, sec=sec)
        Syn[key].tauR       = S_tau_dep
        Syn[key].U          = S_U
        Syn[key].e          = S_e
        
    else:
        sys.stderr.write('\nError: wrong synapse Type (%s). \n\tSynapse not set. Exiting\n' %Type)
        sys.exit()
        
    # create NetStim object
    ns[key]             = h.NetStim()
    ns[key].start       = NS_start
    ns[key].interval    = NS_interval # mean interval between two spikes in ms
    ns[key].noise       = NS_noise
    ns[key].number      = NS_number
    if seed == 'no_seed': pass
    elif seed:  ns[key].seed( seed     )
    else:       ns[key].seed( len(Syn) )
        
    # create NetCon object
    nc[key]             = h.NetCon(ns[key],Syn[key]) #  THIS IS WHERE THE ERROR WAS (Syn[sek] instead of Syn[key])
    nc[key].delay       = NC_delay
    nc[key].weight[0]   = NC_conductance
    nc[key].threshold   = NC_threshold


    
def create_segment_list(cell, dist_groups):
    
    # create segment lists
    segments    = {} 
    
    for d in dist_groups:
        segments[d] = []
          
    # sort segments into list            
    for sec in cell.dendlist:
             
        for seg in sec:
            
            dist = h.distance(seg.x, sec=sec)
            
            if dist < 60:
                segments[0].append(seg)
            elif dist < 120:
                segments[1].append(seg)
            elif dist < 180:
                segments[2].append(seg)
            else:
                segments[3].append(seg)
    
    return segments
    
    

def set_pointers(cell, pointer, mod_list):
    
    for sec in h.allsec():
        for seg in sec:
            for mech in seg:
                
                # if mech in mod_list (skipping car since not used dynamically)
                if mech.name() in mod_list[0:7]:
                    
                    h.setpointer(pointer, 'pka', mech )
    
        
def get_group_and_pattern_index( norm, randActPat ):
    
    
    a = len(randActPat[0])
    b = len(randActPat[1])
    c = len(randActPat[2])
    
    if norm < a:
        group   =   0
        C       =   0
    elif norm < a+b: 
        group   =   1
        C       =   a
    elif norm < a+b+c: 
        group   =   2
        C       =   a+b
    else:
        group   =   3
        C       =   a+b+c
        
    return group, norm-C

    

# alternate version below:
#   set_bg_noise_with_flags
def set_bg_noise(cell,              \
                 cell_type='D1',    \
                 syn_fact=False,    \
                 gabaMod=False,     \
                 fglut=12.0,        \
                 fgaba=4.0,         \
                 dendOnly=0,        \
                 delays=[]          ):
    
    if dendOnly:    compartments = cell.dendlist
    else:           compartments = cell.allseclist
    
    ns      = {}
    nc      = {}
    Syn     = {}
    for s,sec in enumerate(compartments):
        
        # set bg noise----------------------------------
        
        if cell_type == 'D1':
            gbase = 0.3e-3
        else:
            gbase = 0.2e-3
        
        if len(delays) == 0:
            delay = 0
        else:
            delay = delays[s]
            
        # create a glut synapse (glutamate)
        random_synapse(ns, nc, Syn, sec, 0.5,           \
                                NS_interval=1000.0/fglut,    \
                                NC_conductance=gbase,       \
                                NS_start=delay,             \
                                seed=None )
        # create a gaba synapse (Exp2Syn)
        random_synapse(ns, nc, Syn, sec, 0.1,           \
                                Type='gaba',                \
                                NS_interval=1000.0/fgaba,       \
                                NC_conductance=gbase*5,     \
                                NS_start=delay,             \
                                seed=None      )
        
        Syn[sec.name()+'_glut'].ratio = 1.0
        
        if syn_fact:
            Syn[sec.name()+'_glut'].ampa_scale_factor = syn_fact[0]
            Syn[sec.name()+'_glut'].nmda_scale_factor = syn_fact[1]
            
        
        if gabaMod:
            # scale gaba
            nc[sec.name()+'_gaba'].weight[0] = gbase * 3 * gabaMod
        
    
    return Syn, nc, ns

# more intricate copy
def set_bg_noise_with_flags(cell,              \
                 cell_type='D1',    \
                 syn_fact=False,    \
                 gabaMod=False,     \
                 delays=[],         \
                 seedHolder=None,   \
                 skip_compartment={}):
    
    ns      = {}
    nc      = {}
    Syn     = {}
    flag    = { 'soma':{'gaba':True, 'glut':True},
                'axon':{'gaba':True, 'glut':True},
                'dend':{'gaba':True, 'glut':True}}
                
    # update flags based on argument (where and what types of synapses to place)
    for compartment in skip_compartment:
        flag[compartment] = skip_compartment[compartment]
        
    for s,sec in enumerate(cell.allseclist):
        
        # get compartment 
        comp = sec.name().split('[')[0].split('.')[0]
        
        # set bg noise----------------------------------
        
        if cell_type == 'D1':
            gbase = 0.3e-3
        else:
            gbase = 0.2e-3
        
        if len(delays) == 0:
            delay = 0
        else:
            delay = delays[s]
        
        if flag[comp]['glut']:    
            
            # update seed?
            seed=None
            if seedHolder:
                if comp in seedHolder:
                    seed = seedHolder[comp]['glut']
                    seedHolder[comp]['glut'] += 1
            print(seed)    
            # create a glut synapse (glutamate)
            random_synapse(ns, nc, Syn, sec, 0.5,           \
                                    NS_interval=1000.0/12.0,    \
                                    NC_conductance=gbase,       \
                                    NS_start=delay,             \
                                    seed=seed )
            Syn[sec.name()+'_glut'].ratio = 1.0
            if syn_fact:
                Syn[sec.name()+'_glut'].ampa_scale_factor = syn_fact[0]
                Syn[sec.name()+'_glut'].nmda_scale_factor = syn_fact[1]
            
        if flag[comp]['gaba']:
            
            # update seed?
            seed=None
            if seedHolder:
                if comp in seedHolder:
                    seed = seedHolder[comp]['gaba']
                    seedHolder[comp]['gaba'] += 1
                    
            # create a gaba synapse (Exp2Syn)
            random_synapse(ns, nc, Syn, sec, 0.1,           \
                                    Type='gaba',                \
                                    NS_interval=1000.0/3,       \
                                    NC_conductance=gbase*3,     \
                                    NS_start=delay,             \
                                    seed=seed )
            if gabaMod:
                # scale gaba
                nc[sec.name()+'_gaba'].weight[0] = gbase * 3 * gabaMod
    
        
    return Syn, nc, ns


def set_ramping_stimuli(cell,               \
                        random_delays,      \
                        index=None,         \
                        cell_type='D1',     \
                        syn_fact=False,     \
                        gabaMod=False,      \
                        low=800,            \
                        high=1300,          \
                        seed='no_seed'      ):
    
    ns      = {}
    nc      = {}
    Syn     = {}
    name2id = {}
    for s,sec in enumerate(cell.dendlist):
        
        # set bg noise----------------------------------
        
        if cell_type == 'D1':
            gbase = 0.3e-3
        else:
            gbase = 0.2e-3
        
        if index is None:
            if not seed == 'no_seed': np.random.seed(seed=s+1)
            delay = np.random.uniform(low=low, high=high)
        else:
            delay = random_delays[index][s]
        
        name2id[sec.name()] = sec
            
        # create a glut synapse (glutamate)
        random_synapse(ns, nc, Syn, sec, 0.5,           \
                                NS_interval=1000.0/12.0,    \
                                NC_conductance=gbase,       \
                                NS_start=delay,             \
                                seed=seed )
        # create a gaba synapse (Exp2Syn)
        random_synapse(ns, nc, Syn, sec, 0.1,           \
                                Type='gaba',                \
                                NS_interval=1000.0/3,       \
                                NC_conductance=gbase*3,     \
                                NS_start=delay,             \
                                seed=seed      )
        
        Syn[sec.name()+'_glut'].ratio = 1.0
        
        if syn_fact:
            Syn[sec.name()+'_glut'].ampa_scale_factor = syn_fact[0]
            Syn[sec.name()+'_glut'].nmda_scale_factor = syn_fact[1]
            
        
        if gabaMod:
            # scale gaba
            nc[sec.name()+'_gaba'].weight[0] = gbase * 3 * gabaMod
        
    
    return Syn, nc, ns, name2id
    
    
def create_id_mapper():
    '''
    create dict for mapping node index to model.
    For model removal in dspn (71 best)
    '''    
    id_mapper = {}
    counter = 0
    for model_version in range(71):
        if model_version in [43, 8, 66, 54, 31, 32, 57, 45, 30, 25, 68, 67, 19, 21, 53, 6, 60]: continue
        id_mapper[counter] = model_version
        counter += 1
    return id_mapper

def draw_random_factors_DA_dspn( mod_list, modulate_kaf=1 ):
    ''' 
    updated 2019-11-26 from ranges in list -> dict
       this way mod list does not have to be sorted in the exact same order as ranges
       and factors can be drawn for select channels only
    '''
    mod     = {}; ['naf', 'kas', 'kaf', 'kir','cal12','cal13','can']
    if modulate_kaf:
        # set bounds
        b1=0.75; b2=0.85
    else: b1=0.95; b2=1.05
    ranges   = {'naf':  [0.6,0.8],
                'kas':  [0.65,0.85],
                'kaf':  [b1,b2],
                'kir':  [0.85,1.25],
                'cal12':[1.0,2.0],
                'cal13':[1.0,2.0],
                'can':  [0.0,1.0]   }
    for i,channel in enumerate(mod_list):
        mod[channel] = np.random.uniform( ranges[channel][0],ranges[channel][1] )
    return mod

def draw_random_factors_ACh( mod_list ):
    '''
    factors are multiplied onto channel conductance, except for the kaf channel.
    For the kaf channel the gates are shifted into more hyperpolarized domains (0-20 mV)
    '''
    mod     = {}
    ranges  = { 'naf':  [1.0,1.2],
                'kaf':  [0.0,10.0],
                'kir':  [0.8,1.0],
                'cal12':[0.3,0.7],
                'cal13':[0.3,0.7],
                'can':  [0.65,0.85],
                'Im':   [0.0,0.4]   }
    for channel in mod_list:
        mod[channel] = np.random.uniform( ranges[channel][0],ranges[channel][1] )
    return mod

def draw_random_factors_DA_ispn( mod_list ):
    '''
    factors are later multiplied onto channel conductance.
    '''
    mod     = {}
    ranges  = { 'naf':  [0.95,1.05],
                'kaf':  [1.0,1.1],
                'kas':  [1.0,1.1],
                'kir':  [0.8,1.0],
                'cal12':[0.7,0.8],
                'cal13':[0.7,0.8],
                'can':  [0.9,1.0],
                'car':  [0.6,0.8]   }
    for channel in mod_list:
        mod[channel] = np.random.uniform( ranges[channel][0],ranges[channel][1] )
    return mod
       

def set_channel_modulation( org_chan_gbar,          \
                            mod_list,               \
                            mod_fact,               \
                            todict=True,            \
                            modulate_axon=False     ):
    
    '''set modulation of ion channels'''
    
    if todict:
        # transform factors to dict
        factors = {}
        for i,chan in enumerate(mod_list):
            # calc factor
            factors[chan] = mod_fact[i]
    else:
        factors = mod_fact
        
    
    for sec in org_chan_gbar:
        
        # modulate axon?
        if sec.name().find('axon') >= 0:
            if not modulate_axon:
                continue
        
        for seg in org_chan_gbar[sec]:
            
            for mech in org_chan_gbar[sec][seg]: 
                    
                    # get factor from list
                    modulation = factors[ mech.name() ]
                    
                    if mech.name()[0] == 'c':
                        # Ca channels
                        mech.pbar   = org_chan_gbar[sec][seg][mech] * modulation
                        
                    else:
                        # non Ca channels
                        mech.gbar   = org_chan_gbar[sec][seg][mech] * modulation




def make_list_of_gbar(cell, mod_list):
    # makes a list of gbar values in all segments and all mechanims in mod_list.
    # you could skip the sec layer if you like and only save segments, like org_chan_gbar[seg] = ...
    
    org_chan_gbar = {}
    
    for sec in cell.allseclist:
        org_chan_gbar[sec] = {}
        for seg in sec:
            org_chan_gbar[sec][seg] = {}
            for mech in seg:
                
                if mech.name() in mod_list:
                    
                    if mech.name()[0] == 'c':
                        org_chan_gbar[sec][seg][mech] = mech.pbar
                    else:
                        org_chan_gbar[sec][seg][mech] = mech.gbar
    
    return org_chan_gbar



def mixed_stimuli_annealing(activation_pattern, p, s, stim, ncon, rand, spike_list):
    
    # move spikes from donor synapse to acceptors synapse and set the donor activation to []
    
    # spike according to order in list
    pattern     = activation_pattern['patterns' ][p]
    donors      = activation_pattern['donors'   ][p]
    acceptors   = activation_pattern['acceptors'][p]
    steps       = activation_pattern['steps'    ] 
    
    
    if s > 0:
        c       = int(np.sum(steps[:s]))
    else: c     = 0
        
    for i in range(steps[s]):
        
        # find keys to donor and acceptor
        d   =   pattern[ donors[   c+i] ]
        a   =   pattern[ acceptors[c+i] ]
        
        # update spike trains and list
        del stim[a]
        stim[a] = h.VecStim()
        v       = sorted( spike_list[a] + spike_list[d] )
        vec     = h.Vector(v)
        stim[a].play(vec)
        
        spike_list[a] = v
        
        # re-create NetCon object
        ncon[a]            = h.NetCon(stim[a], rand[a])
        ncon[a].delay      = 1
        ncon[a].weight[0]  = (h.synaptic_strength/1000.0)*1e-3 # (uS) = 1.5 nS  Distdep ((1000.0-d2soma)/1000.0)*1e-3
        
        del ncon[d]
        del stim[d]
        del rand[d]
        
        
    
    return spike_list, stim, ncon
        
               



def set_mixed_stimuli(cell, activation_pattern, i, ISI=1, delta=0, syn_fact=False):
    
    rand={}; stim={}; ncon={}; cur={}; spike_list={'name_to_sec':{}}
    
    # spike sequence according to order in list
    pattern     = activation_pattern['patterns' ][i]
    
    interv      = ISI*100
    delay       = np.arange(1000+delta, 1000+delta+interv, ISI)
    
    # loop over the dendritic sections and check set synapse if in list
    for sec in cell.dendlist:
        for section in pattern:
            
            # if section already in spike_list -> continue
            if section in spike_list: continue
            
            if sec.name() == 'dend['+str(int(section))+']':
                
                # get index of activation
                index = np.where(pattern==section)[0].tolist()  # pattern.tolist().index(section)
                
                # create synapse
                rand[section]        =   h.glutamate(0.5, sec=sec) # changed from x -> seg.x
                rand[section].ratio  =   1.0/3.0
                
                if syn_fact:
                    rand[section].ampa_scale_factor = syn_fact[0]
                    rand[section].nmda_scale_factor = syn_fact[1]


                # create VecStim object
                stim[section] = h.VecStim()
                vec = h.Vector(delay[index])
                stim[section].play(vec)
                
                spike_list[section] = delay[index].tolist()     # [delay[index]]

                # create NetCon object
                ncon[section]            = h.NetCon(stim[section], rand[section])
                ncon[section].delay      = 1
                ncon[section].weight[0]  = (h.synaptic_strength/1000.0)*1e-3 # (uS) -> 1.0 nS  Distdep ((1000.0-d2soma)/1000.0)*1e-3 
                
                # map setion name to pointer
                spike_list['name_to_sec'][section] = sec
                
    
    return rand, ncon, stim, spike_list



def set_mixed_stimuli_inSpine(cell, SPINES, pattern, name2sec, 
                                                        ISI=1,  
                                                        syn_fact=False, 
                                                        start_time=1000,
                                                        gmax_nmda=1000,
                                                        ampa2nmda_ratio=1.0/3.0):
    STIM = {}; spike_list={'name_to_sec':{}}
    
    interv      = ISI*100
    delay       = np.arange(start_time, start_time+interv, ISI)
    
    # loop over sections in pattern
    for index,section in enumerate(pattern):
        STIM[section] = {}
        
        # section
        sec     = name2sec[section]  
        
        # attache spine in middle of section
        SPINES[index].connect_spine(sec, x=0.5)
        
        # create synapse, vecstim and netcon objects
        STIM[section]['rand']           =   h.glutamate(0.5, sec=SPINES[index].head) 
        STIM[section]['rand'].ratio     =   ampa2nmda_ratio
        
        if syn_fact:
            STIM[section]['rand'].ampa_scale_factor = syn_fact[0]
            STIM[section]['rand'].nmda_scale_factor = syn_fact[1]
        
        STIM[section]['vstim'] = h.VecStim()
        vec = h.Vector([delay[index]])
        STIM[section]['vstim'].play(vec)
        
        STIM[section]['ncon']            = h.NetCon(STIM[section]['vstim'], STIM[section]['rand'])
        STIM[section]['ncon'].delay      = 0
        STIM[section]['ncon'].weight[0]  = (gmax_nmda/1000.0)*1e-3    # (uS) 
        
        # update spike list
        spike_list[section] = delay[index].tolist()     # [delay[index]]
        spike_list['name_to_sec'][section] = sec
        
    return STIM, spike_list
                
        
                
    
    
        
    
def reset_mixed_stimuli(spike_list, ISI=1, delta=0, syn_fact=False):
    
    rand={}; stim={}; ncon={}
    
    name2sec = spike_list['name_to_sec']
    
    interv      = ISI*100
    delay       = np.arange(1000+delta, 1000+delta+interv, ISI)
    
    # for spike_list entry with len larger than 1:
    i = 0
    for sl in spike_list:
        if len(spike_list[sl]) > 1:
            if sl == 'name_to_sec': continue
            
            sec     = name2sec[sl]
            spikeT  = [delay[i]]
            
            # create synapse
            rand[sl]        =   h.glutamate(0.5, sec=sec) 
            rand[sl].ratio  =   1.0/3.0
            
            if syn_fact:
                rand[sl].ampa_scale_factor = syn_fact[0]
                rand[sl].nmda_scale_factor = syn_fact[1]


            # create VecStim object
            stim[sl] = h.VecStim()
            vec = h.Vector(spikeT)
            stim[sl].play(vec)
            
            spike_list[sl] = spikeT     

            # create NetCon object
            ncon[sl]            = h.NetCon(stim[sl], rand[sl])
            ncon[sl].delay      = 1
            ncon[sl].weight[0]  = (h.synaptic_strength/1000.0)*1e-3 # (uS) -> 1.0 nS  Distdep ((1000.0-d2soma)/1000.0)*1e-3 
            
            i += 1
            
    return rand, ncon, stim, spike_list





def set_clustered_stimuli(  cell,               \
                            section,            \
                            n=20,               \
                            act_time=1000,      \
                            syn_fact=False,     \
                            delta=0,            \
                            ISI=1,              \
                            x=0.5               ):
    
    dend_name   = 'dend[' + str(int(section)) + ']'
    
    for sec in cell.allseclist:
                                                    
        if sec.name() == dend_name:
            
            # calc distance to soma
            d2soma = int(h.distance(x, sec=sec))
            
            # define synapse
            syn         = h.glutamate(x, sec=sec)
            syn.ratio   = 1.0/3.0
            
            if syn_fact:
                syn.ampa_scale_factor = syn_fact[0]
                syn.nmda_scale_factor = syn_fact[1]


            # create NetStim object
            stim            = h.NetStim()
            stim.number     = n
            stim.start      = act_time+delta
            stim.interval   = ISI # mean interval between two spikes in ms (default 1 ms)
            

            # create NetCon object
            ncon             = h.NetCon(stim, syn)
            ncon.delay       = 0
            ncon.weight[0]   = (h.synaptic_strength/1000.0)*1e-3 # (uS). default 1.5 nS
            
            break 
    
    return syn, stim, ncon, d2soma
    


def set_single_E_stim(  section,            \
                        n=1,                \
                        act_time=1000,      \
                        gmax_nmda=1000,     \
                        ratio=1.0/3.0,      \
                        syn_fact=False,     \
                        ISI=1,              \
                        x=0.5               ):
    '''
    sets and returns a excitatory synaptic stimuli in a single section (section).
    
    default values:
    -n = 1              -> number of activations
    -act_time = 1000    -> time of first activation (ms)
    -gmax_nmda = 1000   -> maximal conductance of nmda component (nS; ampa is gmax_nmda * ratio)
    -ratio = 1.0/3.0    -> ratio of ampa/nmda (ampa gmax = gmax_nmda * ratio)
    -syn_fact = False   -> list including factors for ampa and nmda scaling (modulation)
    -ISI = 1            -> interspike interval (ms), averge time between spikes if n > 1
    -x = 0.5            -> location in section, used for chosing activation segment 
    '''
    #TODO: add possibility to use vecStim instead of NetStim
    
    # define synapse
    syn         = h.glutamate(x, sec=section)
    syn.ratio   = ratio
    
    if syn_fact:
        syn.ampa_scale_factor = syn_fact[0]
        syn.nmda_scale_factor = syn_fact[1]


    # create NetStim object
    stim            = h.NetStim()
    stim.number     = n
    stim.start      = act_time
    stim.interval   = ISI # mean interval between two spikes in ms (default 1 ms)
    

    # create NetCon object
    ncon             = h.NetCon(stim, syn)
    ncon.delay       = 0
    ncon.weight[0]   = (gmax_nmda/1000.0)*1e-3 # (uS). default 1.0 nS
    
    return {'syn':syn, 'stim':stim, 'ncon':ncon}


def set_dispersed_stimuli(  cell,               \
                            activationPattern,  \
                            ISI=1,              \
                            syn_fact=False,     \
                            delta=0,            \
                            gabaMod=False       ):
    
    interv  = ISI*10
    
    rand    = {}
    stim    = {}
    ncon    = {}
    vmL     = {}
    cur     = {}
    detect  = {}
    delay   = np.arange(1000+delta, 1000+delta+interv, ISI)
    
    # set random stimuli
    for sec in cell.dendlist:
        
        for seg in sec:  
        
            if seg in activationPattern:
                
                key = sec.name() + '_' + str(seg.x)
                
                # get index of activation
                index = activationPattern.index(seg)
                
                # create synapse
                rand[key]        =   h.glutamate(seg.x, sec=sec) # changed from x -> seg.x
                rand[key].ratio  =   1.0/3.0
                
                if syn_fact:
                    rand[key].ampa_scale_factor = syn_fact[0]
                    rand[key].nmda_scale_factor = syn_fact[1]


                # create NetStim object
                stim[key]            = h.NetStim()
                stim[key].number     = 2
                stim[key].start      = delay[index]
                stim[key].interval   = interv # interval between two spikes in ms (e.g. ISI = 2 -> 1000 / 50 Hz = 20 ms)

                # create NetCon object
                ncon[key]            = h.NetCon(stim[key], rand[key])
                ncon[key].delay      = 1
                ncon[key].weight[0]  = (h.synaptic_strength/1000.0)*1e-3 # (uS) = 1.5 nS  Distdep ((1000.0-d2soma)/1000.0)*1e-3 
                
                #vmL[key] = h.Vector()
                #vmL[key].record(sec(seg.x)._ref_v)
                cur[key] = h.Vector()
                cur[key].record(rand[key]._ref_I)
                
                # set pointer to right target...
                pointer             = rand[key]._ref_i_ampa
    
    return cur, rand, ncon, stim 





def get_fixed_modulation_factors(cell_type, mod_list, modulation_percentage):
    
    
    maximum_modulation = {'D1': {   'naf':0.75,\
                                    'kas':0.8,\
                                    'kaf':0.8,\
                                    'kir':1.2,\
                                    'cal12':1.3,\
                                    'cal13':1.3,\
                                    'can':0.5,\
                                    'car':1.0,\
                                    'ampa':1.2,\
                                    'namda':1.3,\
                                    'gaba':0.8 },   \
                          'D2': {   'naf':1.2,\
                                    'kas':1.1,\
                                    'kaf':1.0,\
                                    'kir':0.9,\
                                    'cal12':0.8,\
                                    'cal13':1.3,\
                                    'can':1.0,\
                                    'car':0.7,\
                                    'ampa':0.8  ,\
                                    'namda':1.0,\
                                    'gaba':0.6 }    }
    
    factors = []
    for i,chan in enumerate(mod_list):
        
        # calc factor
        factors.append( 1 + (maximum_modulation[cell_type][chan]-1) * modulation_percentage / 100 )
    
    return factors
        
        
                            
                            



def analyse_spike_data(spikes, max_time=200, t2s_cutoff=2):
    
    # count spikes and calculate variance of time to first spike and spike ratio
    #
    # spikes = array with size (10,20)
    # max_time = max time after pattern activation for counting of spikes
    # t2s_cutoff = min number of non nan (spikings) over backgrounds
    
    C1              = np.zeros((10))         # mean count after pattern trigger
    C2              = np.zeros((10))         # mean count before pattern trigger
    Ratio           = np.zeros((10))         # mean ratio after pattern trigger
    time2spike      = np.zeros((10))
    
    spiking_traces  = 0
    
    
    # loop over bg
    for b,bg in enumerate(range(10)):
        
        # one row at the time
        spike_train  =   spikes[b,:] 
        
        # index and value of spikes within range (0-max_time)
        I = [ i for i,s in enumerate(spike_train) if s > 0 and s < max_time]
        
        
        # check if spikes
        if len(I) > 0:
            
            spiking_traces += 1
            
            i   = I[0]                          # index of first spike
            s   = spike_train[i]                # first spike (after activation of pattern)
            
        else:
            i   =   0
            s   =   np.nan
            
            
            
        # save number of spikes in range, number of spikes spikes before activation and 
        #      time to first spike after pattern activation
        C1[ b]          =   len(I)
        C2[ b]          =   len([ S[0] for S in enumerate(spike_train) if S[1] < 0])
        time2spike[b]   =   s
        
        
    # calc mean and variance
    Ratio   =  spiking_traces / 10.0 
    
    
    # calc number of non nan values in time2spike
    L       = time2spike.size - np.count_nonzero(np.isnan(time2spike))
    
    if L >= t2s_cutoff:
        Var_t2s  =  np.nanvar(time2spike) 
        mean_t2s =  np.nanmean(time2spike)
    else:
        Var_t2s  =  np.nan
        mean_t2s =  np.nan
                
    
    return [Ratio,          \
            np.mean(C1),    \
            np.mean(C2),    \
            Var_t2s,        \
            mean_t2s        ]




def check_sliding_average(trace, time, reverse=False, split=40, threshold=-40):
    
    # create sliding average and check if larger than threshold (40 mV) over an extended time period
    N       = len(trace)
    n       = int(N/split)
    
    if reverse:
        index = list(range(N-n,0,-n))
    else:
        index = list(range(n,N,n))
    
    x_slide = []
    y_slide = []
    flag    = False
    counter = 0         # keep track of consecutive index over threshold
    
    for s in index:
        e = s+n
        y = np.mean(trace[s:e])
        
        y_slide.append(y)
        x_slide.append(time[s])
        
        if y > threshold:
            
            counter += 1
            
            if counter > 4:
                if reverse:
                    return True
                else:
                    flag = True
        
        else: counter = 0
    
    if reverse:
        return False
    else:
        return [flag, x_slide, y_slide]


# ===== Morphological manipulations     ==================================================

def move_subtree(donor, acceptor):
    '''
    Uses NEURON functions to moves a subtree from one section to another
    
    donor_sec = nrn section object; parent to subtree (all child sections will be moved)
    acceptor_sec = nrn section object where to attache the subtree
    '''
    
    for child_sec in donor.children():     
        h.disconnect(sec=child_sec)
        child_sec.connect(acceptor(1),0)
    

def create_name2sec_dict():
    name2sec = {}
    for sec in h.allsec():
        name = sec.name()
        n = name.split('[')[1].split(']')[0]
        name2sec[name] = sec
        name2sec[n] = sec
    return name2sec


def split_section(cell, sec, n_sec_per_um=1):
    '''
    splits a section into multiple sequentially connected sections.
    
    Retains:
    -total surface area
    -approximate channel distribution
    '''
    
    #import pdb
    #pdb.set_trace()()
    
    # get section values (area etc)
    nseg    = sec.nseg
    L       = sec.L
    area    = 0
    sec_mech = {}
    for i,seg in enumerate(sec):
        area += seg.area()
        sec_mech[i] = {}
        for mech in seg:
            if mech.name() in cell.dendritic_channels:
                if mech.name()[0] == 'c':   sec_mech[i][mech.name()] = mech.pbar
                else:                       sec_mech[i][mech.name()] = mech.gbar
                
    
    # get parent and child sections
    children = []
    for child_sec in sec.children():
        children.append(child_sec)
    parent = sec.parentseg().sec
        
        
    # DISCONNECT section from parent and child sections (starting with child sections).
    for child_sec in children:  
        h.disconnect(sec=child_sec)
    h.disconnect(sec=sec)    
    
    
    # CREATE sequentially new sections
    N = int(np.floor(L/n_sec_per_um))
    if N%2 == 0:
        N += 1  # make sure odd number of sections (and therefore segments)
    
    SEC = []
    for i in range(N):
        sub_sec_name = 'spineSec_%s_%d' % (sec.name(), i)
        ssec = h.Section(name=sub_sec_name)
        SEC.append( ssec )
    
    # SET related values of section to retain overall structure
    diam    = area / (np.pi * L)
    A       = area / N
    for i,ssec in enumerate(SEC):
        ssec.nseg       = 1
        ssec.L          = 1.0 * L / N
        ssec(0.5).diam  = diam
        ssec.insert('pas')
        ssec.Ra         = sec.Ra
        ssec.cm         = sec.cm
        ssec.e_pas      = sec.e_pas
        ssec.g_pas      = sec.g_pas
        value_mapper    = np.floor(nseg * i / N) # values taken from "closest" seg
        for mech in sec_mech[value_mapper]:
            ssec.insert(mech)
            if mech[0] == 'c':  cmd = 'ssec(0.5).pbar_%s = %g' % (mech, sec_mech[value_mapper][mech])
            else:               cmd = 'ssec(0.5).gbar_%s = %g' % (mech, sec_mech[value_mapper][mech])
            exec(cmd)
        ssec.ena        = sec.ena
        ssec.ek         = sec.ek
        ssec.insert("cadyn")
        ssec.insert("caldyn")
    
    # CONNECT new sections
    for i in range(1,N):
        SEC[i].connect(SEC[i-1](1),0)
    # -add child sections
    for child_sec in children:
        child_sec.connect(SEC[-1](1),0)
    
    # -connect new subtree to parent
    SEC[0].connect(parent(1),0)
    
    return SEC
    
    
        

# ===== Functions for vizualisation     ==================================================   


def make_colormap(s, N, CM='magma', custom_map=False):
    # make colormap for viz. of range variable onto morph.
    #
    # s: shape/plotshape object
    # N: number of colors wanted in map
    # CM: python colormap to use. Set to False to use custom map
    # Custom_map: list of list with rgb values, E.g. [ [256,256,256], [0,0,256] ]
    #             the custom map should have length = len(N)
    
    s.colormap(N, 1)
    
    if CM:
        from matplotlib import cm
        cmap    =   cm.get_cmap(CM).colors
        parts   =   len(cmap) / (N-1)
    else:
        cmap    =   custom_map
        parts   =   1
        
    for i in range(N):
        j = parts*i
        print(j)
        r = cmap[j][0]
        b = cmap[j][1]
        g = cmap[j][2]
        s.colormap(i, r,g,b)
        

def get_color_index(value,              \
                    min_max=[-80,40],   \
                    N_colors=13         ):
    # map value of variabel to color index. Used for vizualization of range variable over 
    #       morphology.
    
    # check if value within range 
    if value < min_max[0]:
        return 0
    elif value > min_max[1]:
        return N_colors-1
    else:
        
        # divide the range of values into windows (size based on N_colors)
        n       = (min_max[1]-min_max[0]) / (N_colors-1)
        windows = np.arange(min_max[0],min_max[1]+n/2,n)
        
        # map the value onto a window
        for i in range( N_colors ):
            if np.abs( value - windows[i] ) <= (n/2):
                return i
    
    
    
    

def make_PlotShape( diameter_style=0,   \
                    axis_range=170,     \
                    range_variable="v", \
                    m1=-80, m2=40,      \
                    points=False        ):
    # make plotshape window
    # to use the mechanism spikeMonitor use range_variable="spike_spikeMonitor"
    
    s = h.PlotShape()
    s.show(diameter_style)
    s.size(-axis_range, axis_range, -axis_range, axis_range)
    s.view(-axis_range, -axis_range, axis_range*2,axis_range*2, 500, 500, 1500, 1500)
    s.variable(range_variable)
    s.scale( m1, m2 )
    if points:
        for point in points:
            s.point_mark(2, "O", 10, sec=point)
    s.exec_menu("Shape Plot")
    return s
    

def make_shape( diameter_style=0,   \
                axis_range=170,    \
                range_variable="v", \
                m1=-80, m2=40,      \
                points=False,       \
                zoom=False          ):
    # make shape window
    # TODO update zoom (and pointers?)
    
    s = h.Shape()
    s.show(diameter_style)
    s.size(-axis_range, axis_range, -axis_range, axis_range)
    #s.scale( m1, m2 )
    if zoom:
        # create zoom window.   s.view(mleft, mbottom, mwidth, mheight, sleft, stop, swidth, sheight)
        s.view(-size, -size, size*2,size*2, 500, 500, 1500, 1500)
    if points:
        for point in points:
            s.point_mark(2, "O", 10, sec=point)
    s.exec_menu("Shape Plot")
    return s

        
# ===== USEFULL CODE FROM OTHER SCRIPTS ==================================================  


# how to fit function to data
'''
    # REGRESSION
    
    from scipy.optimize import curve_fit
    
    # sort lists
    I = [(v,i) for i,v in enumerate(D)]
    print(len(I))
    X = [x for _,x in sorted(zip(I,D))]
    Y = np.divide( [x for _,x in sorted(zip(I,R))], norm )
    
    # get regression lines
    popt, pcov = curve_fit(func, X, Y, p0=[1,40,-27])
    
    print(popt)
    
    
    # plot
    ax1.plot(X, Y, '--', color='brown', lw=2)
    ax1.plot(X, func(X, *popt), 'k', lw=6) '''
    



# how to show structure morphology
'''
    h.toplogy()'''


# how to record current
'''
    # record current
    channels = ['naf', 'kas', 'kaf', 'kir', 'kdr', 'sk', 'bk', 'cal12', 'cal13', 'can', 'car', 'cat32', 'cat33']
    rec      = {}
    SEG      = {}
    
    for sec in h.allsec():
        Type = sec.name().split('[')[0]
        if Type in SEG:
            continue
        
        for seg in sec:

            if Type in SEG:
                break
            else:
                SEG[Type] = True

            for mech in seg:

                if mech.name() in channels:

                    # build record vector
                    name = mech.name() + '_' + Type
                    rec[name] = h.Vector()
                    if name[0] == 'k':
                        rec[name].record(mech._ref_gk)
                    elif name[0] == 'n':
                        rec[name].record(mech._ref_gna)
                    elif name[2] == 'l':
                        rec[name].record(mech._ref_gcal)
                    else:
                        rec[name].record(mech._ref_gca)'''
    

# example code how to partly block channels
'''
    # configure simulation to record from both calcium pools.
    # the concentration is here summed, instead of averaged. 
    # This doesn't matter for the validation fig, since relative concentration is reported.
    # For Fig 5B, where concentration is reported, this is fixed when plotting
    # (dividing by 2 since the volume of both pools are the same).
    # -> see the plot_Ca_updated function in plot_functions.
    for i,sec in enumerate(h.allsec()):
        
        for j,seg in enumerate(sec):
            
            #print( sec.name(), sec.name().split('[')[1].split(']')[0])            
            
            # uncomment here if testing kaf blocking effect on bAP
            if block == 'kas':
                
                # blocking kas (100%)
                seg.kas.gbar   = 0
                
            elif block == 'naf': # and sec.name().find('the right section') >= 0:
        
                # block naf (100%) for all segments 30-50 mu from soma 
                dist = h.distance(seg.x, sec=sec)
                if dist > 10 and dist < 60:
                    seg.naf.gbar   = 0
            
            
            if sec.name().find('axon') < 0: # don't record in axon
                
                if block == 'kaf':
                
                    # blocking kaf (50% ?), kaf is not distributed to the axon
                    block_fraction = 0.999
                    gbar           = seg.kaf.gbar
                    seg.kaf.gbar   = (1 - block_fraction) * gbar
                
                sName = sec.name().split('[')[0]
                
                
                # N, P/Q, R Ca pool
                cmd = 'ca_%s%s_%s = h.Vector()' % (sName, str(i), str(j))
                exec(cmd)
                cmd = 'ca_%s%s_%s.record(seg._ref_cai)' % (sName, str(i), str(j))
                exec(cmd)   
                
                # the L-type Ca
                cmd = 'cal_%s%s_%s = h.Vector()' % (sName, str(i), str(j))
                exec(cmd)
                cmd = 'cal_%s%s_%s.record(seg._ref_cali)' % (sName, str(i), str(j))
                exec(cmd)  
                
    # solver------------------------------------------------------------------------------            
    cvode = h.CVode()
    
    h.finitialize(cell.v_init)
    
    # run simulation
    while h.t < tstop:
                
        h.fadvance()
    
            
    # save output ------------------------------------------------------------------------
    
    res         = {}

    if cell_type == 'D1':
        distances   = np.arange(40,200, 10)
    else:
        distances   = np.arange(49,200, 10)


    # for all sec
    for i,sec in enumerate(h.allsec()):
        # skip axon
        if sec.name().find('axon') < 0:
            # for all seg in sec
            for j,seg in enumerate(sec):

                dist    =    int( np.round( h.distance(seg.x) ) ) 
                
                # if dist is out of range
                if dist <= (distances[0]  - 5) or \
                   dist >= (distances[-1] + 5):
                        
                        continue
                        
                # calc index (k) and get distance key (d) 
                x       =    (dist - distances[0])/10.0
                k       =    int(round(x))
                d       =    distances[k]
                
                # save to result dict
                if d not in res:
                    res[d] = []

                vName   =   'ca_%s%s_%s'  %  ( sName, str(i), str(j)  )
                v2Name  =   'cal_%s%s_%s' %  ( sName, str(i), str(j)  )
                cmd     =   'V = np.add(%s, %s)' % (vName, v2Name) # this is were concentrations are summed (see above ~line 140)
                exec(cmd)
                dCa     =   max(V[3300:-1]) - V[3300]
                res[d].append(dCa)
    
    # calc mean and return list as function of distance
    y       = []
    norm    = np.mean(res[distances[0]])
    
    for d in distances:
        y.append( np.divide(np.mean(res[d]), norm) )

    # return vector             
    return y    '''   
    
