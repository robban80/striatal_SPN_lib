

import numpy as np
import pickle, glob
from scipy.signal import butter, filtfilt, freqz
import matplotlib.pyplot as plt


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y 

def extract_features_local(t,v, threshold=-20):
    
    # get index of values larger than threshold (-) ....---....---....--- 
    index_above_threshold = [x[0] for x in enumerate(v) if x[1] > threshold]
    
    if len(index_above_threshold) < 1:
        return []
    
    # get first index of each spike sequence
    index = []
    for i in range(1,len(index_above_threshold)):
        
        if not index_above_threshold[i] == index_above_threshold[i-1]+1:
            index.append( index_above_threshold[i] )
            #print(i, index_above_threshold[i], index_above_threshold[i-1])
    
    
    # loop over all spike sequences and extract max values
    extreme_values  = { 'v' :{  'max':np.zeros(len(index)+1,dtype=int), 
                                'min':np.zeros(len(index)+1,dtype=int)},
                        'i' :{  'max':np.zeros(len(index)+1,dtype=int),
                                'min':np.zeros(len(index)+1,dtype=int)}}
    start_index     = index_above_threshold[0]
    for ii,i in enumerate(index):
        end_index   = i
        spike_seq   = v[start_index:end_index]
        extreme_values['v']['max'][ii] = max(spike_seq)
        extreme_values['v']['min'][ii] = min(spike_seq)
        extreme_values['i']['max'][ii] = np.argmax(spike_seq) + start_index
        extreme_values['i']['min'][ii] = np.argmin(spike_seq) + start_index
        start_index = i
    extreme_values['v']['max'][-1] = max(v[start_index:])
    extreme_values['v']['min'][-1] = min(v[start_index:])  
    extreme_values['i']['max'][-1] = np.argmax(v[start_index:start_index+400]) + start_index
    extreme_values['i']['min'][-1] = np.argmin(v[start_index:start_index+400]) + start_index 
    return extreme_values
     

    

def get_spike_times(t,curve, return_index=False):
    
    # differentiate
    dtrace = np.diff(curve)
        
    # take sign
    sign = np.sign(dtrace)
    
    # diff sign change
    dsign = np.diff(sign)
    
    # find pos/neg values
    index = np.add(np.where(dsign < 0.0), 1)   # np.argwhere(dsign > 0.0) + 1
    
    # get time points
    spike_times = {'t':t[index].tolist()[0], 'y':curve[index]}
    
    if return_index:
        return spike_times, index
    else:
        return spike_times


def get_current_direction_from_soma(vsoma, vaxon, vdend_dict):
    '''
    Calculates the direction of current flow between soma and neurites by comparing
        the somatic membrane potential against the initial section of each neurite
    -vsoma and vaxon are list with membrane potential over time 
    -vdend is a dict containing lists of all dendritic potentials
    * returns a dict with differences over time 
    '''
    
    # initial sections number of each dendritic stem
    initial_sections = [0, 19, 30, 31, 42, 49, 54, 55]
    
    # loop over inital sections
    DF = {'axon':np.subtract(vsoma,vaxon)} 
    for seck in vdend_dict:
        secID = int(seck.split('[')[1].split(']')[0])
        if secID in initial_sections:
            # calc diff and add to list
            i = initial_sections.index(secID)
            DF['dend%d'%(i)] = np.subtract(vsoma,vdend_dict[seck])
            
    return DF
    
def vary_mgblock(ax=None, mg=1.0, alpha=-0.062, beta=3.57):
    
    if not ax: fig,ax = plt.subplots(1,1)
    
    v = np.arange(-90,50)
    mgblock =  1.0 / (1 + mg * np.exp(alpha * v) / beta )
    ax.plot(mgblock, v, label=str(alpha))
    

def plot_mgblock(ax, c_mg=0, c_2nd=0, p_1stDer=1, p_fbtw=1, minV=-80, maxV=40,
                 alpha=-0.062, beta=3.57, return_peak=False):
    
    v = np.arange(minV,maxV)
    mgblock =  1.0 / (1 + 1.0 * np.exp(alpha * v) / beta )
    diff = np.diff(mgblock)
    ac = np.diff(diff)
    if p_fbtw:ax.fill_between([0,1], [-42,-42], [-47,-47], color='k', alpha=0.5, lw=0)
    if p_1stDer: ax.plot(np.divide(diff,max(diff)), v[:-1])
    if c_2nd: ax.plot(np.divide(ac,max(ac)), v[1:-1], color=c_2nd)
    else:     ax.plot(np.divide(ac,max(ac)), v[1:-1])
    if c_mg:  ax.plot(mgblock, v, color=c_mg)
    else:     ax.plot(mgblock, v)
    
    maxind = np.argmax(ac)
    ax.plot( [0,0.1], [v[1:-1][maxind],v[1:-1][maxind]], ls='-', color=c_2nd )
    
    ax.set_xlim([1.05,0])
    
    if return_peak:
        return v[1:-1][maxind]
    
    
    
def plot_window_current(ax, color='k', minV=-80, maxV=40, return_peak=True):
    
    v           = np.arange(minV,maxV)
    mVhalf      = -25.0 
    hVhalf      = -62.0 
    mSlope      =  -9.2 
    hSlope      =   6.0     
    minf        = 1 / (1 + np.exp( (v-mVhalf) / mSlope ) )
    hinf        = 1 / (1 + np.exp( (v-hVhalf) / hSlope ) )
    
    window      = np.multiply( np.power(minf,3), hinf )
    wnorm       = np.divide(window, max(window))
    ax.plot( wnorm, v, ls='--', color=color )
    maxind = np.argmax(wnorm)
    print(maxind,wnorm[maxind], v[maxind])
    ax.plot( [0,0.1], [v[maxind],v[maxind]], ls='-', color=color )
    
    if return_peak:
        return v[maxind]



def best_fit_parameter_distribution(array, parameters, ax=None):
    '''
    plot distribution range of best fit solutions (normalized to parameter range)
    '''
    
    ax = ax if ax is not None else plt.gca()
    
    colors = ['k', 'orange']
    alphas = [0.3, 1.0]
    
    # parameter ranges
    param_range      = [    [-0.5,0.5], \
                            [0.8,1.0], \
                            [10.0,60.0],   \
                            [1.0,30.0],  \
                            [-5,5], \
                            [-5,5], \
                            [0.1, 0.58], \
                            [-0.5,0.5], \
                            [-0.5,0.5], \
                            [0.0,0.9], \
                            [1.0,130.0],   \
                            [-70.0,-3.0],  \
                            [-0.5,0.5],    \
                            [-5.0,60.0],  \
                            [1.0,70.0], \
                            [-0.5,0.5], \
                            [-0.5,0.5],    \
                            [-9.0,-6.0], \
                            [1.0,130.0],   \
                            [-70.0,-3.0], \
                            [-9.0,-6.0], \
                            [1.0,130.0],   \
                            [-70.0,-3.0], \
                            [-7.0,-5.0], \
                            [0.8,1.0], \
                            [10.0,60.0],   \
                            [1.0,30.0] ]
    
    #fd,ad   =   plt.subplots(1,1, figsize=(16,8))
    
    x       =   list(range( len(parameters)-2))
    
    # loop rows in array
    for index in range(array.shape[0]):
        
        color = colors[int(array[index,1])]
        alpha = alphas[int(array[index,1])]
        
        # create list
        y = []
        for j in x:
            val     =   array[index,j+2]
            A       =   param_range[j][0]
            B       =   param_range[j][1]
            factor  =   (val-A) / (B-A)
            y.append(   factor  )
            
            if val > B or val < A:
                print( j, parameters[j+2], A, val, B ) 
                print()
            
        # plot
        ax.plot(x, y, '-o', ms=20, color=color, alpha=alpha)
             
    ax.set_xticks(x)
    ax.set_xticklabels(parameters[2:], fontsize=20, rotation=90)



def check_sliding_average(data, return_cs_index=False, threshold=-40):
    '''
    low pass filter and...
    checks if sliding mean is over -40 mV for at least x locations in the trace (should perhaps be consequtive).
    If so it is classified as complex. 
    Not very robust but works for traces that spikes sparsely.
    
    len trace should be 500 ms and dt 25um.
    
    splits trace in 40 pieces and checks mean of each piece sequentially.
    '''
    
    # low pass filter parameters        
    order   = 2
    fs      = 15.0      # sample rate, Hz
    cutoff  = 0.01      # desired cutoff frequency of the filter, Hz
    # filter local potential
    trace = butter_lowpass_filter(data, cutoff, fs, order)
    
    N = len(trace)
    n = int(N/40)
    
    over_index = []
    
    #print len(trace), n
    if not return_cs_index:
        for s in range(N-n,0,-n):
            e = s+n
            if np.mean(trace[s:e]) > -40:
                over_index.append(s)
                if len(over_index) > 3:
                    return True
    else:
        for i in range(0,N,n):
            e = i+n
            if np.mean(trace[i:e]) > threshold:
                over_index.append(i)
                if len(over_index) > 2:
                    return True, over_index, trace
    if return_cs_index:           
        return False, over_index, trace
    else: return False

def check_sliding_average_lowSampleData(data, return_cs_index=False, threshold=-40):
    '''
    low pass filter and...
    checks if sliding mean is over -40 mV for at least 4 locations in the trace (should perhaps be consequtive).
    If so it is classified as complex. 
    Not very robust but works for traces that spikes sparsely.
    
    len trace should be 500 ms and dt 25um.
    
    splits trace in 40 pieces and checks mean of each piece sequentially.
    '''
    
    # low pass filter parameters        
    
    order   = 2
    fs      = 2.0      # sample rate, Hz
    cutoff  = 0.004      # desired cutoff frequency of the filter, Hz
    # filter local potential
    trace = butter_lowpass_filter(data, cutoff, fs, order)
    
    N = len(trace)
    n = int(N/120)
    
    over_index = []
    
    # check trace from end to beginning
    for s in range(N-n,0,-n):
        e = s+n
        if np.mean(trace[s:e]) > threshold:
            over_index.append(s)
            if len(over_index) > 1:
                if return_cs_index:
                    return True, over_index, trace
                else:
                    return True
        else:
            over_index = []
    if return_cs_index:           
        return False, over_index, trace
    else: 
        return False

def check_sliding_average_lowSampleData_from0(data, return_cs_index=False, threshold=-40):
    '''
    low pass filter and...
    checks if sliding mean is over -40 mV for at least 4 locations in the trace (should perhaps be consequtive).
    If so it is classified as complex. 
    Not very robust but works for traces that spikes sparsely.
    
    len trace should be 500 ms and dt 25um.
    
    splits trace in 40 pieces and checks mean of each piece sequentially.
    '''
    
    # low pass filter parameters        
    
    order   = 2
    fs      = 2.0      # sample rate, Hz
    cutoff  = 0.004      # desired cutoff frequency of the filter, Hz
    # filter local potential
    trace = butter_lowpass_filter(data, cutoff, fs, order)
    
    N = len(trace)
    n = int(N/120)
    
    over_index = []
    
    # check trace from end to beginning
    for s in range(0,N,n):
        e = s+n
        if np.mean(trace[s:e]) > threshold:
            over_index.append(s)
            if len(over_index) > 1:
                if return_cs_index:
                    return True, over_index, trace
                else:
                    return True
        else:
            over_index = []
    if return_cs_index:           
        return False, over_index, trace
    else: 
        return False

def check_tract_data(time, trace):

    # find index for t > 0 (clustered activation time) 
    index = next(x[0] for x in enumerate(time) if x[1] >= 0)   
    
    # ---- check if complex spikes following clustered activation
    if check_sliding_average(trace[index:]):
        cs = 1
    else: cs = 0
    
    # ---- average vm before clustered synaptic activation
    
    # clip potential spikes (to not affect average)
    clipped_vm = [-50 if x > -50 else x for x in trace[:index]]
    
    mean_vm = np.mean( clipped_vm ) 
    
    
    return mean_vm, cs
    

def loop_over_pickled_files(files):

    # loop over pickle files
    for f in files:
        
        # get data from file name
        if f.find('Axon') >= 0:
            axon    = True
            c       = count['axon']
        else: 
            axon = False
            c       = count['all']
        
        if c > 18: continue
        
        run     = f.split('run')[  1].split('_')[0]
        model   = f.split('model')[1].split('_')[0]
        
        # load file and extract information
        with open(f, 'rb') as handle:
            data = pickle.load(handle)
        
        # -holds complex spikes?
        vrest, cs = check_tract_data(data['time'], data['Vm']['soma'])
        
        # -get naf parameter data
        if axon:    random_variables = params_axon_only
        else:       random_variables = params_all
        
        
        shiftm = random_variables[c]['variables']['naf_shift'][0]
        shifth = random_variables[c]['variables']['naf_shift'][1]
        taum   = random_variables[c]['variables']['naf_shift'][2]
        tauh   = random_variables[c]['variables']['naf_shift'][3]
        
        chan_params = [vrest, cs]
        for chan in parameters:
            for p in range(len(random_variables[c]['variables'][chan])):
                chan_params.append( random_variables[c]['variables'][chan][p] )
        
        # sort into pandas frame
        if axon:
            array_axon[count['axon'],:] = chan_params
            count['axon'] += 1
            i = 0
        else:
            array_all[ count[ 'all'],:] = chan_params
            count[ 'all'] += 1
            i = 1
        
        # plot
        ax[i][cs].plot( data['time'], data['Vm']['soma'], 'k' )


def extract_extreme_values(i,f, regions=['soma']):	
        
    try:
        with open(f, 'rb') as handle:
            result_dict = pickle.load(handle)
    except:
        return None
    
    t  = np.array([x for x in  result_dict['time'] if x >= 0])
    
    if i == 0:
        all_data    = {'file':f, 'time':t}
    else: all_data  = {'file':f}
    
    for r,region in enumerate(regions):
        
        # classify
        if region == 'soma':
            v = np.array([result_dict['Vm'][region][x[0]] for x in  enumerate(result_dict['time']) if x[1] >= 0])
            mean_vm, cs = check_tract_data(t, v)
        elif region == 'axon':
            v = np.array([result_dict['Vm'][region][x[0]] for x in  enumerate(result_dict['time']) if x[1] >= 0])
        
        if region == 'dendrites':
            all_data[region] = {}
            if region not in result_dict['Vm']: continue
            for dend in result_dict['Vm'][region]:
                
                v = np.array([result_dict['Vm'][region][dend][x[0]] for x in  enumerate(result_dict['time']) if x[1] >= 0])
                # get extreme values (peaks and dipp)
                extreme_values = extract_features_local(t,v,threshold=-25)
                # add to dict
                all_data[region][dend] = {'cs': cs, 'extreme_values': extreme_values, 'vm':v}
        else:   
            
            # get extreme values (peaks and dipp)
            extreme_values = extract_features_local(t,v,threshold=-25)
            # add to dict
            all_data[region] = {'cs': cs, 'extreme_values': extreme_values, 'vm':v}

    return all_data



def extract_extreme_values2(result_dict, region='soma'):
    '''
    checks if complex spike and if so extracts extreme values before cs. 
    Differences compared to first version are
        - only traces classified as complex are handled.
        - traces are trucated at cs (if cs), i.e. no data after cs are kept.
        - only checks one cellular region
    '''
    
    
    # remove all points before stimuli
    t       = np.array([x for x in  result_dict['time'] if x >= 0])
    trace   = np.array([result_dict['Vm'][region][x[0]] for x in enumerate(result_dict['time']) if x[1] >= 0])
    # check if complex
    cs, ind, lpft = check_sliding_average(trace, return_cs_index=True)
    
    if not cs: 
        # get index of min (and max) peaks 
        maxi = next(x[0] for x in enumerate(t) if x[1] >= 120)
        extreme_values = extract_features_local(t,trace[:maxi],threshold=-5)
        return None, ind, trace, lpft, extreme_values
    else:
        # truncate trace at complex spike
        maxi = np.argmax( lpft )
        # get index of min (and max) peaks 
        extreme_values = extract_features_local(t,trace[:maxi],threshold=-5)
        return cs, ind, trace, lpft, extreme_values
        # synchronize traces on peaks before spikes

def extract_extreme_general(trace, time, threshold=0.5):
    '''get extreme values by looking at changes in derivative
    
    filter by taking peaks with amplitude of minimum half min of threshold*min(trace)'''
    
    # diff
    d = np.diff(trace)
    
    # sign 
    s = np.sign(d)
    
    # min
    min_peak_index = (np.diff(s) > 0).nonzero()[0]+1
    
    # max
    max_peak_index = (np.diff(s) < 0).nonzero()[0]+1
    
    # filter to only use ap peaks within 120 ms
    peaks = np.array( [i for i in min_peak_index if trace[i] < threshold*min(trace) and time[i] < 120] )
    
    # first max after last min (last max before cs)
    lmcs = next( i for i in max_peak_index if i > peaks[-1] )
    
    return {'min':peaks, 'max':lmcs}
    



def extract_cs_dipps(trace, time, csi=None):
    '''uses changes in derivative
    filter by taking peaks with amplitude of minimum half min of threshold*min(trace)'''
    
    if not csi: 
        maxt = 500
        mint = 0
    else:
        maxt = time[csi]
        mint = maxt-150
        
    # diff
    d = np.diff(trace)
    
    # sign 
    s = np.sign(d)
    
    # min
    min_peak_index = (np.diff(s) > 0).nonzero()[0]+1
    
    # filter to only use ap peaks within 120 ms
    peaks = np.array( [i for i in min_peak_index    if trace[i] < 0.55*min(trace) 
                                                    and time[i] < maxt 
                                                    and time[i] > mint
                                                    and trace[i] > -50 ])
    
    # max
    max_peak_index = (np.diff(s) < 0).nonzero()[0]+1
    '''
    plt.plot( time, trace)
    plt.plot(np.array(time)[peaks], np.array(trace)[peaks], '-ok')
    plt.plot(np.array(time)[csi], np.array(trace)[csi], 'or')
    plt.show()'''
    
    # first max after last min (last max before cs)
    lmcs = next( i for i in max_peak_index if i > peaks[-1] )
    
    return {'min':peaks, 'max':lmcs}



def extract_cs_dur(trace, time, start_index):
    '''returns dur of cs (symetric from start_index)'''
    
    print(start_index)
    
    # start time and voltage
    tstart = time[start_index]
    vstart = trace[start_index]
    
    # end time
    Vm = trace[start_index:]
    Tm = time[start_index:]
    
    tend = next( Tm[i] for i,v in enumerate(Vm) if v < vstart )
    
    # used to create example of cs extraction
    '''
    plt.figure()
    plt.plot([tstart, tend], [vstart, vstart], lw=2, zorder=2)
    plt.plot(Tm, Vm, 'k', lw=3)
    plt.plot(time[start_index], trace[start_index], 'o', zorder=4)
    
    plt.axis('off')
    plt.savefig('./Figures/example_csDurExtraction.png', dpi=300, transparent=True)
    plt.show()'''
    
    return tend-tstart
    
   
def plot_Ina(a, tShift, region, result_dict, color):
    j       = ['axon','soma','dend'].index(region)
    
    ina = np.array([result_dict['na'][region]['I'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0])
    a[2,j].plot( tShift, ina, color=color )
    if region =='axon':
        # peak values of naf in the AIS
        ev = extract_extreme_general(ina, tShift, threshold=0.5)
        a[2,j].plot(  tShift[ev['min']], ina[ev['min']], 'o', color=color )
        a[2,j].plot(  tShift[ev['max']], ina[ev['max']], 'o', color=color )
        t2shift = tShift + 100 - tShift[ ev['min'][-1] ]
        a[2,3].plot( t2shift, ina, color=color )
        a[2,3].plot(  t2shift[ev['min']], ina[ev['min']], 'o', color=color )
        a[2,3].plot(  t2shift[ev['max']], ina[ev['max']], 'o', color=color )
    else: t2shift = None        
    a[3,j].plot(  tShift, 
                [result_dict['na'][region]['m'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0],
                color=color   )
    a[4,j].plot(  tShift, 
                [result_dict['na'][region]['h'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0],
                color=color   )    
    return t2shift


                
def plot_nmda(a, tShift, t2Shift, result_dict, color):
    #if 'nmda' in result_dict:
    if not 'nmda' in result_dict: return
    if 'I' in result_dict['nmda']:
        for nmda in list(result_dict['nmda']['I'].values()):
            n = [nmda[ii] for ii,x in  enumerate(result_dict['time']) if x >= 0]
            a[0,3].plot( tShift,  n, color=color )   
            a[1,3].plot( t2Shift, n, color=color )
    else:
        for nmda in list(result_dict['nmda'].values()):
            n = [nmda[ii] for ii,x in  enumerate(result_dict['time']) if x >= 0]
            a[0,3].plot( tShift,  n, color=color )   
            a[1,3].plot( t2Shift, n, color=color )


def calc_mean_vm_dend(result_dict, time):
    sum_dend = np.array(time, dtype='float64')
    for sec,dend in result_dict['Vm']['dendrites'].items():
        adend = np.array([dend[ii] for ii,x in  enumerate(result_dict['time']) if x >= 0])
        sum_dend += adend
    return sum_dend / len(result_dict['Vm']['dendrites'])    



def create_name2secDict():
    ''' 
    OBS, not tested since realized not needed...
    creates and returns a map from dendritic section name to neuron section object
    '''
    import neuron as nrn
    import MSN_builder as build 
    cell = build.MSN(  params='../D1-MSN_wip/params_dMSN.json',                  
                       morphology='../D1-MSN_wip/WT-dMSN_P270-20_1.02_SGA1-m24.swc'  )
    name2sec = {}
    for dend in cell.dendlist:
        name2sec[dend.name()] = dend
    return name2sec
        
    

def plot_dendritic_vm_vs_dist(result_dict, index, ax, c):
    
    with open('../D1-MSN_wip/Libraries/map_sec_dist_to_soma.pkl', 'rb') as f:
        dist_mapper = pickle.load(f)
    
    index -= int(0/0.025)
    
    # collect
    dist = np.zeros(len(dist_mapper))
    vm   = np.zeros(len(dist_mapper))   
    i    = 0 
    for sec,dend in result_dict['Vm']['dendrites'].items():
        dist[i] = dist_mapper[ sec.split('[')[1].split(']')[0] ] 
        # TODO: shift index in call to function -> remove loop
        vm[i]   = [dend[ii] for ii,x in  enumerate(result_dict['time']) if x >= 0][index]
        i += 1
    
    # sort and plot
    index = np.argsort(dist)
    ax.plot(dist[index], vm[index], color=c)
    
    
        
        
        

                   
def load_and_plot(i,a):
    colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf']
    files = glob.glob('Pickled_recordings/*.pkl')
    
    region='soma'
    
    # open file    
    try:
        with open(files[i], 'rb') as handle:
            result_dict = pickle.load(handle)
    except:
        print('-failed to load', i)
        return False,False,False,False
    
    time = np.array([x for x in  result_dict['time'] if x >= 0])
       
    try:
        cs, ind, trace, lpft, extr = extract_extreme_values2(result_dict, region=region)
    except:
        print('-failed to extract', i)
        return False,False,False,False
    '''
    a[0,0].plot(time, trace, c=colors[i%8])
    a[0,0].plot(time, lpft, c=colors[i%8]) 
    a[0,0].plot(time[ind], lpft[ind], '-o', c=colors[i%8], lw=3)
    a[0,0].plot(time[extr['i']['min']], trace[extr['i']['min']], '-o', c=colors[i%8])
    a[0,0].plot(time[extr['i']['max']], trace[extr['i']['max']], '-o', c=colors[i%8])
    '''
    if 'i' in extr:
        len_spikes = len(extr['i']['max'])
        if len_spikes > 1:
            Imax = extr['i']['max']
            isi = [time[Imax[ii+1]]-time[iii] for ii,iii in enumerate(Imax[:-1])]
            print (isi)
        else: isi=False
        if 'dendrites' in result_dict['Vm']:
            if cs: c = 'r'
            else:  c = 'k'
            plot_dendritic_vm_vs_dist(  result_dict, 
                                        extr['i']['min'][-1],
                                        a[3,3],
                                        c)
    else:
        len_spikes = 0
        isi = False
    if cs:  
        print('--', i)                                                                       
        tShift = time + 100 - time[ extr['i']['min'][-1] ] 
        master_index = extr['i']['min'][-1]                         
        a[0,0].plot(    tShift, 
                        [result_dict['Vm']['axon'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0], 
                        c=colors[i%8])
        a[0,1].plot(    tShift, 
                        [result_dict['Vm']['soma'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0], 
                        c=colors[i%8]) 
        a[0,1].plot(tShift[master_index], trace[master_index], '-o', c=colors[i%8])
        a[0,1].plot(tShift[extr['i']['max']], trace[extr['i']['max']], '-o', c=colors[i%8])
        # sodium
        t2Shift = plot_Ina(a, tShift, 'axon', result_dict, colors[i%8])
        blaj    = plot_Ina(a, tShift, 'soma', result_dict, colors[i%8])
        plot_nmda(a, tShift, t2Shift, result_dict, colors[i%8])
        # k+ axon
        '''
        a[2,1].plot(  tShift, 
                    [result_dict['k']['M'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0],
                    color=colors[i%8]   )
        a[3,1].plot(  tShift, 
                    [result_dict['k']['kas'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0],
                    color=colors[i%8]   )'''
        # dend
        if 'dendrites' in result_dict['Vm']:
            mean_dend = calc_mean_vm_dend(result_dict, time)
            a[0,2].plot(tShift, mean_dend, c=colors[i%8])
        return True, len_spikes, isi, extr
    else:
        a[1,1].plot(time, trace, c=colors[i%8])
        a[1,1].plot(time, lpft,  c=colors[i%8])
        a[1,1].plot(time[ind], lpft[ind], 'o', c=colors[i%8])
        a[1,0].plot(    time, 
                        [result_dict['Vm']['axon'][ii] for ii,x in  enumerate(result_dict['time']) if x >= 0], 
                        c=colors[i%8])
        # dend
        if 'dendrites' in result_dict['Vm']:
            mean_dend = calc_mean_vm_dend(result_dict, time)
            a[1,2].plot(time, mean_dend, c=colors[i%8])
        return False, len_spikes, isi, extr

def synchronize_complex_spikes(n):
    plt.close('all')
    f,a = plt.subplots(5,4, figsize=(24,24))
    #columns axon|soma|dend
    # rows   vm-cs, vs-normal, Ina, mna, hna
    b = 1
    all_isi_befor_cs = {}
    CS_data = {'nspikes':[], 'lastISI':[], 'allISI':{}, 'extrSpk':{} }
    for i in range((b-1)*100,777):    # b*100
        #if not i in [0, 3, 10, 13, 16, 18, 24, 27, 28, 37, 38, 42, 44, 47, 51, 55, 56, 58, 61, 63, 67, 72, 74, 76, 78, 83, 90, 91, 93, 99, 109, 111, 114, 120, 126, 132, 140, 142, 144, 148, 158, 159, 162, 163, 166, 168, 196, 198, 202, 205, 212, 221, 222, 223, 227, 228, 236, 238, 247, 250, 258, 262, 266, 267, 268, 269, 272, 274, 275, 279, 280, 282, 287, 305, 315, 318, 332, 334, 335, 340, 341, 347, 350, 356, 357, 358, 359, 362, 363, 367, 377, 379, 381, 389, 393, 394, 395, 396, 397, 400, 403, 404, 412, 418, 426, 433, 435, 436, 440, 447, 460, 466, 478, 480, 485, 493, 500, 501, 509, 511, 513, 516, 519, 520, 522, 523, 532, 538, 539, 540, 554, 556, 559, 560, 561, 562, 568, 571, 587, 593, 596, 601, 602, 607, 612, 614, 615, 627, 634, 644, 645, 647, 655, 659, 660, 661, 662, 670, 677, 679, 692, 693, 703, 704, 713, 714, 718, 719, 720, 723, 732, 745, 749, 756, 758, 762, 766, 770, 776]: continue
        cs, ls, isi, extr = load_and_plot(i,a)
        if cs:
            CS_data['nspikes'].append(ls)
            CS_data['extrSpk'][i] = extr
            if isi:
                CS_data['lastISI'].append(isi[-1])
                CS_data['allISI'][i] = isi
            
    regions = ['axon','soma','dend']
    ylabels = ['vm cs','vm normal','I na', 'm na', 'h na']
    
    # plot windows setup    
    for j in range(3):
        a[0,j].set_title(regions[j])
        for i,ax in enumerate(a[:,j]):
            ax.set_xlim([0,200])
            if j == 0:
                a[i,j].set_ylabel(ylabels[i])  
            if i < 2:
                a[i,j].set_ylim([-80,50])    
    f.savefig('sort_and_synchronize_complex_spikes_all.png') 
    #plt.close('all')
    plt.figure()
    plt.hist(CS_data['nspikes'], bins=5, range=(0.5,5.5))
    plt.xlim(0,6)
    plt.title('spikes before cs')
    plt.savefig('sort_and_synch_hist_spikes_before_cs_all.png')
    plt.figure()
    plt.hist(CS_data['lastISI'], bins=40)
    plt.title('last isi before cs')
    plt.savefig('sort_and_synch_hist_last_isi_before_cs_all.png')
    #plt.close('all')  
    save_lib = False    # this saves cs_data to file, uncomment if rerun (overwrites stored data!)
    if save_lib:  
        with open('lib_CSdata.pkl', 'wb') as handle:
            pickle.dump(    CS_data, 
                            handle, 
                            protocol=pickle.HIGHEST_PROTOCOL)
    plt.show()


def show_result_dict_structure():
    '''does nothing. only here to show key structure in result dict
    will trigger errors if run since none of the key-value pairs are compleate
    '''
    
    all_traces = {  'time':t,
                    'Vm': { 'soma':v, 'axon':va, 'dendrites':ALLDEND},
                    'na': { 'axon':{ 'I':nf,  'm':nm,  'h':nh,  'mtau':mtau,  'htau':htau},
                            'soma':{ 'I':nfs, 'm':nms, 'h':nhs, 'mtau':mtaus, 'htau':htaus}},
                    'k':  { 'M':m, 'kas':ks},
                    'bg': { 'soma_i':bg_i_s, 'soma_e':bg_e_s, 'axon_i':bg_i_a, 'axon_e':bg_e_a},
                    'gaba':{'spikes':spikes, 'sum':gaba_current, 'bg':GABA},
                    'I_syn':il,
                    'gbg':G_BG,
                    'nmda':NMDA,
                    'ampa':AMPA}

# ----------------------------------------------------------------------------------------
# correlate spike features

def hinton(matrix, max_weight=None, ax=None, col=None, row=None):
    """Draw Hinton diagram for visualizing a weight matrix.
    from:
    http://python-for-multivariate-analysis.readthedocs.io/a_little_book_of_python_for_multivariate_analysis.html
    """
    ax = ax if ax is not None else plt.gca()

    if not max_weight:
        max_weight = 2**np.ceil(np.log(np.abs(matrix).max())/np.log(2))

    ax.patch.set_facecolor('lightgray')
    ax.set_aspect('equal', 'box')
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    
    if col and row:
        
        w       = matrix[col][row]
        color   = 'orange' if w > 0 else 'black'
        print( w, color )
        size    = np.sqrt(np.abs(w))
        rect = plt.Rectangle([0.5 - size / 2, 0.5 - size / 2], size, size,
                                 facecolor=color, edgecolor=color)
        ax.add_patch(rect)
        
    else:
        for (x, y), w in np.ndenumerate(matrix):
            color   = 'orange' if w > 0 else 'black'
            size = np.sqrt(np.abs(w))
            rect = plt.Rectangle([x - size / 2, y - size / 2], size, size,
                                 facecolor=color, edgecolor=color)
            ax.add_patch(rect)
    
    
        nticks = matrix.shape[0]
        #ax.xaxis.tick_top()
        ax.set_xticks(list(range(nticks)))
        ax.set_xticklabels(list(matrix.columns), rotation=45, fontsize=30)
        ax.set_yticks(list(range(nticks)))
        ax.set_yticklabels(matrix.columns, fontsize=30)
        ax.grid(False)

        ax.autoscale_view()
        ax.invert_yaxis()
        

def correlate_cs_spike_features():
    import pandas as pd
    
    # import library from file
    with open('lib_CSdata.pkl', 'rb') as handle:
        CS_data = pickle.load(handle)
            
    # extract features
    # CS_data = {'nspikes':[], 'lastISI':[], 'allISI':{}, 'extrSpk':{} }
    
    columns     = [ 'nspikes',  
                    'peak0', 'peak1', 'peak2', 
                    'dipp0', 'dipp1', 'dipp2',
                    'isi0', 'isi1', 'isi2',
                    't_cs']
    nfeatures   = len(columns)
    array       = np.full( [len(CS_data['extrSpk']),nfeatures], np.nan )
    for i,key in enumerate(CS_data['extrSpk']):
        if key in CS_data['allISI']: nisi = len(CS_data['allISI'][key]); 
        else:                        nisi = 0
        array[i,0]  = nisi + 1
        array[i,10] = CS_data['extrSpk'][key]['i']['min'][-1]*0.025 # dt = 0.025 -> time (ms)
        for j in range(1,4):
            if j > nisi:
                array[i,j]      = CS_data['extrSpk'][key]['v']['max'][-j]
                array[i,j+3]    = CS_data['extrSpk'][key]['v']['min'][-j] 
                break
            array[i,j]      = CS_data['extrSpk'][key]['v']['max'][-j]
            array[i,j+3]    = CS_data['extrSpk'][key]['v']['min'][-j]
            array[i,j+6]    = CS_data['allISI' ][key            ][-j]
             
    # import to pandas frame   
    df = pd.DataFrame(array, columns=columns)
        
    # run correlation between different features last spike amp vs #spikes
    corr = df.corr()
    f,ax = plt.subplots(1,1)
    hinton(corr, ax=ax)
    plt.show()
    
