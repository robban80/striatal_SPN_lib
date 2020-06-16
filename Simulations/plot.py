
# plot things for ejn paper

import sys, json, glob, pickle,codecs
#sys.path.insert(0, '../../')
#import common_functions     as use
import numpy                as np
import pandas               as pd
import functions4analysis   as f4a

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# curve_fit function
from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a * np.exp( (x-b) / c )


colors      = {'ctrl':'k','ACh':'#fdc086','DA':'#7fc97f','ACh+DA':'#beaed4'}  
maxT        = 500
minT        = 0
path2traces = 'Results/'



def run():
    
    #plot_all_traces()
    #find_example_traces_ispn()
    fi_curve()                                      # plots panel 6A1 in Lindroos and Hellgren Kotaleski (2020?) 
    bap()                                           # 6A2
    '''plot_example_traces_ispn()                      # 6B1
    ispn_plot_proportions(load=1)                   # 6B2
    ispn_plot_proportions_nafNMDA(load=1)           # 6C1-2
    ispn_plot_cs_factors(load=1)                    # 6C3
    #ispn_plot_proportions_nafNMDA_DA(load=1)'''
 
 
 # --------------------------------------------------------------------------------------------------------------
 
def unserialize_from_json(name, mod=1):
    '''
    Reverse serialization of json into ndarray. To run on json files created using above
    https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable
    '''
    obj_text = codecs.open(name, 'r', encoding='utf-8').read()
    loaded_data = json.loads(obj_text)
    loaded_data['ctrl'] = np.array(loaded_data['ctrl'])
    if mod: loaded_data['mod' ] = np.array(loaded_data['mod' ])
    return loaded_data
      

def fi_curve():
    fig,ax = plt.subplots(1,1, figsize=(5,4) )
    
    color   = ['#969696','#525252']
    color2  = ['#cb181d','#fb6a4a']
    rheob   = {'exp':{1:[],2:[]}, 'mod':{1:[],2:[]}}
    
    # model data --------------------------------------------------------------------
    for i,mtype in enumerate(['dspn','ispn']):
        spikes  = unserialize_from_json('../Validation_data/{}_extracted_spikes_fi_ctrlOnly.json'.format(mtype), mod=0)
        nmodels = [71,34][i]
        
        with open('../Libraries/D{}_{}bestFit_updRheob.pkl'.format(i+1,nmodels), 'rb') as f:
            model_sets = pickle.load(f, encoding="latin1")    
        for model in range(nmodels):
            if mtype == 'dspn' and \
                model in [43, 8, 66, 54, 31, 32, 57, 45, 30, 25, 68, 67, 19, 21, 53, 6, 60]: continue
            rheobase    = model_sets[model]['rheobase']  
            I           = [rheobase-1] + [i+rheobase for i in spikes['I']]
            spike_train = [0] + [s for s in spikes['ctrl'][model,:]]
            ax.plot(I, spike_train, color=color[i],  lw=2, alpha=1)
            
            rheob['mod'][i+1].append(rheobase)
    
    # experimental data ---------------------------------------------------------------
    z = [0, 1e6, 0]
    for d in range(2):
        exp = glob.glob('../Validation_data/Planert2013-D{}-FI-trace*'.format(d+1))
        for i,f in enumerate(exp):
            [x_i,y_i] = np.loadtxt(f, unpack=True)
            ax.plot(x_i, y_i, ls='--', color=color2[d], lw=2, alpha=1, zorder=z[d]+i)
            
            rheob['exp'][d+1].append(x_i[0])
    
    # mean and std rheobase ----------------------------------------------------------
    y1 = -10
    y2 = -5
    ax.plot([np.mean(rheob['exp'][1])-np.std(rheob['exp'][1]), np.mean(rheob['exp'][1])], [y1,y1], '--', c=color2[0], ms=10)
    ax.plot([np.mean(rheob['exp'][2])-np.std(rheob['exp'][2]), np.mean(rheob['exp'][2])], [y2,y2], '--', c=color2[1], ms=10)
    ax.plot([np.mean(rheob['mod'][1])+np.std(rheob['mod'][1]), np.mean(rheob['mod'][1])], [y1,y1], '-', c=color[0], ms=10)
    ax.plot([np.mean(rheob['mod'][2])+np.std(rheob['mod'][2]), np.mean(rheob['mod'][2])], [y2,y2], '-', c=color[1], ms=10)
    
    ax.plot([np.mean(rheob['exp'][1])], [y1], '|', c=color2[0], ms=10)
    ax.plot([np.mean(rheob['exp'][2])], [y2], '|', c=color2[1], ms=10)
    ax.plot([np.mean(rheob['mod'][1])], [y1], '|', c=color[0], ms=10)
    ax.plot([np.mean(rheob['mod'][2])], [y2], '|', c=color[1], ms=10)
    
    ax.plot([np.mean(rheob['exp'][1])], [y1], 'o', c=color2[0], ms=5)
    ax.plot([np.mean(rheob['exp'][2])], [y2], 'o', c=color2[1], ms=5)
    ax.plot([np.mean(rheob['mod'][1])], [y1], 'o', c=color[0], ms=5)
    ax.plot([np.mean(rheob['mod'][2])], [y2], 'o', c=color[1], ms=5)
            
    # frame -----------------------------------------------------------------------------
    ybase=0
    ax.plot([ybase,ybase],[0,25],       'k',lw=3)
    ax.plot([ybase,ybase+10],[25,25],   'k',lw=3)
    ax.plot([ybase,ybase+10],[0,0],     'k',lw=3)
    for i in range(100,601,250):
        ax.plot([i,i],[-0.5,0.5],       'k',lw=3)
    ax.set_xlim([-100,800])
    ax.axis('off')
    fig.savefig('./Figures/validation_FI_freq2.png', transparent=True, dpi=300)


def bap():
    fig, ax     = plt.subplots(figsize=(3.75,5))
    
    with open('../Validation_data/ispn_extracted_bap.json', 'r') as f:
        ispn_res = json.load(f)
    dspn_file = '../Validation_data/dspn_extracted_bap.json'
    with open(dspn_file, 'r') as f:
        dspn_res     = json.load(f)
    
    RES = { 'dspn':{'data':dspn_res, 'min':35, 'max':40}, 
            'ispn':{'data':ispn_res, 'min':45, 'max':50}
            }
    
    color       = ['#969696','#525252']
    color2      = ['#cb181d','#fb6a4a']
    maxdist     = 210
    
    # model data ------------------------------------------------------------------
    for i,mtype in enumerate(['dspn','ispn']):
        distList    = RES[mtype]['data']['dist']
        m1          = RES[mtype]['min']
        m2          = RES[mtype]['max']
        
        # normalize to random (first) compartment within 35-40 um somatic dist
        index       = np.argsort( distList )
        norm_ind,D  = next(i for i in enumerate(distList) if i[1] < m2 and i[1] > m1)
        dist        = [distList[i]  for i in index \
                                    if distList[i] >= D and distList[i] < maxdist]
        for mv in range(71):
            if mtype == 'dspn' and \
                mv in [43, 8, 66, 54, 31, 32, 57, 45, 30, 25, 68, 67, 19, 21, 53, 6, 60]: continue
            elif mtype == 'ispn' and \
                mv > 33: break
            
            caList     = RES[mtype]['data'][str(mv)]
            Csum       = [caList[i]/caList[norm_ind] for i in index \
                                                if distList[i] >= D and distList[i] < maxdist]
            # get regression lines
            popt, pcov = curve_fit( func, 
                                    dist, 
                                    Csum, 
                                    p0=[1,40,-15])
            ax.plot(dist, Csum, 'o', ms=6, c=color[i], mew=1, mec='w', alpha=1)
            ax.plot(dist, func(dist, *popt), color=color[i], lw=1, zorder=10000+mv)
    
    
    # experimental data (extraced from published article) --------------------------------
    for d in range(2):
        [x1,y1]  = np.loadtxt('../Validation_data/bAP-DayEtAl2006-D{}.csv'.format(d+1), unpack=True) 
        ax.plot(x1,y1, color=color2[d], lw=3, ls='--', zorder=1e6)
    
    # figure setup ----------------------------------------------------
    lw=4
    ax.plot([0,5],[0,0],         'k',lw=lw)
    ax.plot([0,5],[1,1],         'k',lw=lw)
    ax.plot([0,0],[0,1],         'k',lw=lw)
    ax.set_xlim([-10,250])
    ax.axis('off')
    fig.savefig('Figures/validation_bap.png', transparent=True, dpi=300)
    
    plt.close('all')   
    
        
def plot_all_traces():
    '''loops trough all files defined by search string and plots included voltage traces
       
       used to get an overview of data
       '''
    
    alpha       = 1.0
    tremove     = 800
    
    for i in range(1,34): 
        
        for c in ['ctrl','ACh','DA','ACh+DA']:
            fstring = 'Results/inVivo_ramping_{}_*_model{}.json'.format(c,i)
            files = glob.glob(fstring)
            for f in files: 
        
                with open(f, 'r') as handle:
                    data = json.load(handle)
                
                fig,ax = plt.subplots(1,1)
                time = [t-1000 for t in data['time'] if t > tremove]
                
                for j in data:
                    if j == 'time': continue
                    for bg in data[str(i)]:
                        if bg == 'par': continue
                        
                        trace = [data[j][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                        ax.plot(time, trace, c=colors[c], alpha=alpha)
                        
                ax.set_title(f)
                plt.show()


def find_example_traces_ispn():
    '''searches input files for complex spikes.
            also plots and prints the traces it finds
       
       used in scanning for cs
       '''
    #ax[0].set_title(mid)
    alpha = 0.5
    
    count = {2:0, 1:0, 0:0}
    
    flag = 0
    
    condition = 'RedNafIncNMDA'
    sim_tag = 'ACh'
    
    files = glob.glob('Results/inVivo_ramping{}_{}_*_model*.json'.format(condition,sim_tag))
    print(len(files))
    fig,ax = plt.subplots(1,3, figsize=(12,4))
    for f in files:
        #f = 'Results/inVivo_ramping_{}_model{}.json'.format(c,mid)
        try:
            with open(f, 'r') as handle:
                data = json.load(handle)
            #print(f)
        except:
            continue
        
        #print(f)
        tremove = 1000
        time = data['time'] #[t-1000 for t in data['time'] if t > tremove]
        # da fid: ach fid: condition: bg id
        for i in data:
            if i == 'time': continue
            for bg in data[str(i)]:
                #print(data[i].keys())
                if bg == 'par': continue
                trace = data[i][str(bg)] #[data[str(i)]['0'][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                
                if max(trace) < -10:
                    #ax[0].plot(time, trace, alpha=alpha)
                    count[0] += 1
                elif f4a.check_sliding_average_lowSampleData(trace, threshold=-37):
                    # check for complex spikes.
                    ax[2].plot(time, trace, alpha=alpha)
                    #plt.plot(time, trace, alpha=alpha)
                    #plt.show()
                    s = "'mod':{}, 'bg':{}".format(i,bg)
                    print('{'+s+'},')
                    count[2] += 1
                    flag=1
                else:
                    #ax[1].plot(time, trace, alpha=alpha)
                    count[1] += 1
    
    print(sim_tag, condition, count)   
             
    for i in range(3):
        ax[i].set_xlabel(count[i])
    if True:
        plt.show()
    else:
        plt.close('all')


def ispn_plot_cs_factors(load=0):
    
    if load:
        with open('Results/cs_factors_ispn.json', 'r') as handle:
            res = json.load(handle)
    else:
        res = {}
        for condition,c in zip(['','RedNaf', 'RedNafIncNMDA', 'IncNMDA'], ['range', 'naf', 'nafNMDA', 'NMDA']):
            count = {2:0, 1:0, 0:0}
            files = glob.glob('Results/inVivo_ramping{}_ACh_*_model*.json'.format(condition))
            print(len(files))
            
            res[c] = []
            
            for f in files:
                
                try:
                    with open(f, 'r') as handle:
                        data = json.load(handle)
                except:
                    print(f)
                
                tremove = 900
                time    = [t-1000 for t in data['time'] if t > tremove]
            
                for j in data:
                    if j == 'time': continue
                    for bg in data[str(j)]:
                        if bg == 'par': continue
                        
                        trace = [data[j][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                        
                        # count spikes and complex spikes
                        if max(trace) < -10:
                            count[0] += 1
                        elif f4a.check_sliding_average_lowSampleData(trace, threshold=-37):
                            count[2] += 1
                            # get factors
                            res[c].append(data[j]['par'])
                            break
                        else:
                            count[1] += 1
        
        # save files
        with open('Results/cs_factors_ispn.json', 'w') as handle:
            json.dump(res, handle, indent=4)
    
    # handle res
    
    for c in res:
        NMDA = []
        NAF  = []
        f1,a1 = plt.subplots(1,1, figsize=(4,4))
        f2,a2 = plt.subplots(1,1, figsize=(4,4))
        for par in res[c]:
            naf  = par['factors']['ach']['intr']['naf']
            nmda = par['factors']['ach']['syn']['NMDA']
            NAF.append(naf)
            NMDA.append(nmda)
            
        a1.hist(NAF,  bins=20, range=(0.7,1.2), color='k' )
        a2.hist(NMDA, bins=20, range=(1.0,1.6), color='k' )
        
        for a,x,x1 in zip([a1,a2],[[0.7,1.2],[1.0,1.6]],[1.2,1.05]):
            a.axis('off')
            a.plot(x,[0,0], 'k', lw=4)
            y = a.axis()
            a.plot([1,1],[0,y[3]], 'lightgrey', lw=2)
            a.fill_between([1,x1],[y[3],y[3]],[0,0], color='lightgrey')
        
        if not c == 'range':
            f1.savefig('Figures/ispn_cs_factors_{}_naf.png'.format(c),  dpi=300, transparent=True)
            f2.savefig('Figures/ispn_cs_factors_{}_nmda.png'.format(c), dpi=300, transparent=True) 
            



def ispn_plot_proportions(load=0):
    '''plot fig 6B2 in Lindroos and Hellgren Kotaleski (2020?)
       
       proportion of spiking traces (and proportion with cs---not in manuscript since no CS is these data)
       values on axes are added in post processing 
       '''
    if load:
        with open('Results/cs_count_.json', 'r') as handle:
            res = json.load(handle)
    else:
        res = {}
        for sim_tag in ['ACh','ACh+DA','ctrl','DA',]:
            count = {2:0, 1:0, 0:0}
            files = glob.glob('Results/inVivo_ramping_{}_*_model*.json'.format(sim_tag))
            print(len(files))
            
            for f in files:
                
                try:
                    with open(f, 'r') as handle:
                        data = json.load(handle)
                except:
                    print(f)
                
                tremove = 900
                time    = [t-1000 for t in data['time'] if t > tremove]
            
                for j in data:
                    if j == 'time': continue
                    for bg in data[str(j)]:
                        if bg == 'par': continue
                        
                        trace = [data[j][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                        
                        # count spikes and complex spikes
                        if max(trace) < -10:
                            count[0] += 1
                        elif f4a.check_sliding_average_lowSampleData(trace, threshold=-37):
                            count[2] += 1
                        else:
                            count[1] += 1
        
            print(sim_tag, count)
            res[sim_tag] = count   
    
        # save files
        with open('Results/cs_count_.json', 'w') as handle:
            json.dump(res, handle, indent=4)
        
    # calc percentage spiking (and CS)
    CS = []
    NS = []
    f1,a1 = plt.subplots(1,1, figsize=(4,2))
    f2,a2 = plt.subplots(1,1, figsize=(4,2))
    for i,cond in enumerate(['ctrl','ACh','DA','ACh+DA']):
        
        ns = res[cond]['1']
        cs = res[cond]['2']
        tot = res[cond]['0'] + res[cond]['1'] + res[cond]['2']
        
        for st,a,n in zip([ns,cs],[a1,a2],['reg','cs']):
            height = 100*st/tot
            print(cond, n, height)
            n1,b1,p1 = a.hist([0,1,2,3], bins=4, color=colors[cond])
            for j,p in enumerate(p1):
                if j==i:
                    p.set_height(height)
                else:
                    p.set_height(0)
    
    for a in [a1,a2]:
        a.plot([0,3],[0,0],'k',lw=2)
        a.plot([0,3],[50,50],'lightgrey',lw=1)
        a.plot([0,3],[100,100],'lightgrey',lw=1)
        a.set_ylim([-2,102])
        a.axis('off')
    f1.savefig('Figures/ispn_proportion_spiking_range.png', dpi=300, transparent=True)
    f2.savefig('Figures/ispn_proportion_cs_range.png', dpi=300, transparent=True)
    plt.show()


def ispn_plot_proportions_nafNMDA(load=0):
    '''plot fig 6C1-2 in Lindroos and Hellgren Kotaleski (2020?)
       
       this data set includes reduced sodium (uniformly) and increased NMDA -> giving CS
            run under ACh condition
            
       proportion of spiking traces (C1) and proportion with cs (C2)
       
       values on axes are added in post processing 
       '''
    fname = 'Results/cs_count_nafNMDA.json'
    
    if load:
        with open(fname, 'r') as handle:
            res = json.load(handle)
    else:
        res = {}
        for condition,c in zip(['','RedNaf', 'RedNafIncNMDA', 'IncNMDA'], ['range', 'naf', 'nafNMDA', 'NMDA']):
            count = {2:0, 1:0, 0:0}
            files = glob.glob('Results/inVivo_ramping{}_DA_*_model*.json'.format(condition))
            print(len(files))
            
            for f in files:
                
                try:
                    with open(f, 'r') as handle:
                        data = json.load(handle)
                except:
                    print(f)
                
                tremove = 900
                time    = [t-1000 for t in data['time'] if t > tremove]
            
                for j in data:
                    if j == 'time': continue
                    for bg in data[str(j)]:
                        if bg == 'par': continue
                        
                        trace = [data[j][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                        
                        # count spikes and complex spikes
                        if max(trace) < -10:
                            count[0] += 1
                        elif f4a.check_sliding_average_lowSampleData(trace, threshold=-37):
                            count[2] += 1
                        else:
                            count[1] += 1
        
            print(c, count)
            res[c] = count   
    
        # save files
        with open(fname, 'w') as handle:
            json.dump(res, handle, indent=4)
        
    # calc percentage spiking (and CS)
    CS = []
    NS = []
    f1,a1 = plt.subplots(1,1, figsize=(4,2))
    f2,a2 = plt.subplots(1,1, figsize=(4,2))
    for i,cond in enumerate(['range', 'naf', 'NMDA', 'nafNMDA']):
        
        ns = res[cond]['1']
        cs = res[cond]['2']
        tot = res[cond]['0'] + res[cond]['1'] + res[cond]['2']
        
        for st,a,n,w in zip([ns,cs],[a1,a2],['reg','cs'],[100,1000]):
            height = w*st/max(tot,1)    # avoid dev by zero for empty set...
            print(cond, n, height)
            n1,b1,p1 = a.hist([0,1,2,3], bins=4, color='grey', edgecolor='w', linewidth=1.2)
            for j,p in enumerate(p1):
                if j==i:
                    p.set_height(height)
                else:
                    p.set_height(0)
    
    for a,w in zip([a1,a2],[100,10]):
        a.plot([0,3],[0,0],'k',lw=2)
        a.plot([0,3],[w/2,w/2],'lightgrey',lw=1)
        a.plot([0,3],[w,w],'lightgrey',lw=1)
        a.set_ylim([-2,w+2])
        a.axis('off')
    f1.savefig('Figures/ispn_proportion_spiking_nafNMDA.png', dpi=300, transparent=True)
    f2.savefig('Figures/ispn_proportion_cs_nafNMDA.png',      dpi=300, transparent=True)
    plt.show()

def ispn_plot_proportions_nafNMDA_DA(load=0):
    '''same as ispn_plot_proportions_nafNMDA, but run under DA conditions (uniform modulation)
            not included in manuscript (fewer complex spikes)
       '''
    #fname = 'Results/cs_count_nafNMDA.json'
    fname = 'Results/cs_count_nafNMDA2_DA.json'
    
    fig,ax = plt.subplots(1,2)
    
    if load:
        with open(fname, 'r') as handle:
            res = json.load(handle)
    else:
        res = {}
        for i,c in enumerate(['DA','ACh+DA']):
            count = {'2':0, '1':0, '0':0}
            files = glob.glob('Results/inVivo_rampingDARedNafIncNMDA_{}_*_model*.json'.format(c))
            print(len(files))
            
            for f in files:
                
                try:
                    with open(f, 'r') as handle:
                        data = json.load(handle)
                except:
                    print(f)
                
                tremove = 900
                time    = [t-1000 for t in data['time'] if t > tremove]
            
                for j in data:
                    if j == 'time': continue
                    for bg in data[str(j)]:
                        if bg == 'par': continue
                        
                        trace = [data[j][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                        
                        # count spikes and complex spikes
                        if max(trace) < -10:
                            count['0'] += 1
                        elif f4a.check_sliding_average_lowSampleData(trace, threshold=-37):
                            count['2'] += 1
                            ax[i].plot(time,trace)
                        else:
                            count['1'] += 1
        
            print(c, count)
            res[c] = count   
    
        # save files
        with open(fname, 'w') as handle:
            json.dump(res, handle, indent=4)
        
    # calc percentage spiking (and CS)
    CS = []
    NS = []
    f1,a1 = plt.subplots(1,1, figsize=(4,2))
    f2,a2 = plt.subplots(1,1, figsize=(4,2))
    for i,cond in enumerate(['DA', 'ACh+DA']):
        print(res[cond].keys())
        ns = res[cond]['1']
        cs = res[cond]['2']
        tot = res[cond]['0'] + res[cond]['1'] + res[cond]['2']
        
        for st,a,n,w in zip([ns,cs],[a1,a2],['reg','cs'],[100,1000]):
            height = w*st/max(tot,1)    # avoid dev by zero for empty set...
            print(cond, n, height)
            n1,b1,p1 = a.hist([0,1,2,3], bins=4, color='grey', edgecolor='w', linewidth=1.2)
            for j,p in enumerate(p1):
                if j==i:
                    p.set_height(height)
                else:
                    p.set_height(0)
    
    for a,w in zip([a1,a2],[100,10]):
        a.plot([0,3],[0,0],'k',lw=2)
        a.plot([0,3],[w/2,w/2],'lightgrey',lw=1)
        a.plot([0,3],[w,w],'lightgrey',lw=1)
        a.set_ylim([-2,w+2])
        a.axis('off')
    ax[0].axis('off')
    fig.savefig('Figures/ispn_cs_DAnafNMDA.png', dpi=300, transparent=True)
    f1.savefig('Figures/ispn_proportion_spiking_DA_nafNMDA.png', dpi=300, transparent=True)
    f2.savefig('Figures/ispn_proportion_cs_DA_nafNMDA.png',      dpi=300, transparent=True)
    plt.show()


def plot_example_traces_ispn(ax=None):
    '''plots fig 6B1 in Lindroos and Hellgren Kotaleski (2020?)
       
       10 first traces under each condition are plotted
       setup: ispn bombarded with ramping synaptic input under concurrent neuromodulation
       '''
       
    if not ax: fig,ax = plt.subplots(1,1)
    
    mid         = 1
    tremove     = 900
    ntraces     = 10
    alpha       = 1
    
    for c in ['ctrl','DA','ACh','ACh+DA']:
        count   = 0
        fstring = 'Results/inVivo_ramping_{}_*_model{}.json'.format(c,mid)
        files   =  glob.glob(fstring)
        for f in files:
            with open(f, 'r') as handle:
                data = json.load(handle)
            
            time = [t-1000 for t in data['time'] if t > tremove]
            for j in data:
                if j == 'time': continue
                for bg in data[str(j)]:
                    if bg == 'par': continue
                    
                    trace = [data[j][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                    ax.plot(time, trace, c=colors[c], alpha=alpha, lw=2, zorder=np.random.randint(1000))
                    
                    count += 1
                    if count >= ntraces:
                        flag = 1
                        break
                break
            break
        
    
    lw=3
    ax.plot([minT-50,maxT], [-70, -70], '--k', lw=2, zorder=0) 
    
    xbase=0; ybase=-40; xlen=-100; ylen=-20
    ax.plot([xbase,xbase+xlen], [ybase,ybase], 'k', lw=lw)
    ax.plot([xbase,xbase], [ybase+ylen,ybase], 'k', lw=lw)
    ax.axis('off')
    fig.savefig('Figures/ispn_example_traces_in_range.png', dpi=300, transparent=True)
    plt.show()

   

run()

