
# plot things for ejn paper

import sys, json, glob, pickle
#sys.path.insert(0, '../../')
#import common_functions     as use
import numpy                as np
import pandas               as pd
import functions4analysis   as f4a


from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#plt.style.use('ggplot')
 
colors      = {'ctrl':'k','ACh':'#fdc086','DA':'#7fc97f','ACh+DA':'#beaed4'}  
csingle     = {'dend':'grey', 'soma':'k', 'aprox':'m', 'adist':'r'}
conditions  = ['ctrl','ACh','DA','ACh+DA']
ontop       = [3,2,1,0]
maxT        = 500
minT        = 0
nbg         = 5 
nmodels     = 54


def run():
    
    find_example_traces_ispn()
    
    # main panel 
    #main_panel_fig_cs_intro()
    
    # detailed analysis of 1 trace
    #detailed_analysis_single_trace()
    
    # inhibition
    #inhibition_fig()
    
    # complex spikes
    #cs_factors()
    
        
        
    


     
    

# -------------------------------------------------------------------------------------
# 1st level


def detailed_analysis_single_trace():
    
    ctrl = 0
    mdl  = 2   # 2, 5, 10 
    
    with open('Results/single_inVivo_ramping_{}_model{}.json'.format(['ACh+DA','DA'][ctrl],mdl), 'r') as handle:
        data = json.load(handle)
    
    if ctrl:    save_tag = 'model{}_ctrl'.format(mdl)
    else:       save_tag = 'model{}'.format(mdl)
    
    
    
    if False:
        # vm
        fig, ax = get_ax()
        plot_single_vm_soma_axon(data, ax=ax)
        fig.savefig('Figures/cs_example_vm-axon_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
        fig, ax = get_ax()
        plot_single_vm_soma_dend(data, ax=ax)
        fig.savefig('Figures/cs_example_vm-dend_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
    
    if False:
        # synaptic currents
        fig, ax = get_ax()
        ax[1].axis('off')
        plot_single_glut(data, ax=ax[0])
        fig.savefig('Figures/cs_example_synaptic_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
        fig,ax   = plt.subplots(2,1, figsize=(4,4),gridspec_kw = {'hspace':0.8})
        correlate_nmdaDend_with_vmSoma(data,ax=ax)
        fig.savefig('Figures/correlation_vm-nmda_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
        plt.close('all')
    
    if True:
        # hist showing cs initiation voltage
        plot_cs_dipps(load=True)
        plot_cs_dur(load=True)
        plot_multiple_nmdas_spikes_cs_dur()
        plt.show()
    
    if False:
        # intrinsic currents (sodium)
        fig, ax = get_ax()
        ax[1].axis('off')
        plot_single_sodium(data, ax=ax[0])
        fig.savefig('Figures/cs_example_Isoma_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
        #plt.show()
        plt.close('all')
    
    if False:
        # Ca
        fig, ax = get_ax()
        ax[1].axis('off')
        plot_single_ca(data, ax=ax[0])
        fig.savefig('Figures/cs_example_Ca_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
        #plt.close('all')
        #plt.show()
    
    if False:
        # diff vm between compartment, showing soma/dendrites ending CS
        fig     = plt.figure(figsize=(6,4))
        r       = 6
        grid    = plt.GridSpec(r, r, hspace=1.0, wspace=0)
        a1      = fig.add_subplot(grid[:3,:-1])
        a2      = fig.add_subplot(grid[3:,:-1])
        a3      = fig.add_subplot(grid[:,-1])
        a3.axis('off')
        ax = [a1,a2]
        plot_single_vm_diff(data, ax=ax)
        fig.savefig('Figures/cs_example_diffVm_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
        
        # enlarged sodium showing enlarged current
        fig     = plt.figure(figsize=(6,4))
        r       = 6
        grid    = plt.GridSpec(r, r, hspace=1.0, wspace=0)
        a1      = fig.add_subplot(grid[:3,:-1])
        a2      = fig.add_subplot(grid[3:,:-1])
        a3      = fig.add_subplot(grid[:,-1])
        a3.axis('off')
        ax = [a1,a2]
        plot_single_vm_sodium_zoomed(data, ax=ax)
        fig.savefig('Figures/cs_example_Isoma_zoom_ramping_{}.png'.format(save_tag), dpi=300, transparent=True)
    
    plt.show()


def cs_factors():
    ''' analyse extracted data in pandas data frame (df)
    '''
    
    #fig,ax  = plt.subplots(1,1,figsize=(5,6))
    # plot 
    #   -proportion of cs over condition
    #   -proportion of spikes over condition
    #   -proportion of cs over model (ACh+DA)
    #       + correlation: rheobase, nafAUC, kafAUC
    
    uniform = 0
    if uniform:
        path    =   'Results/UniformDAmod/'
        tag     =   '_uniform'
        subs    =   'UniformMod'
    else:
        path    =   'Results/'
        tag     =   ''
        subs    =   ''
        
    #df  = extract_data_into_frame(load=1, path=path)
    #df_proportion_condition(df)
    #df_proportion_model(df, tag='_uniform') #, ax=ax)
    #example_traces_factors()
    cs_correlation_with_naf(load=0, subset=subs) # subset: ''=no axon; 'UniformMod'=uniform
    
    '''
    fig, ax = plt.subplots(1,1,figsize=(6,3))
    df_proportion_cs(df, ax=ax)
    fig.savefig('Figures/cs_proportion{}_ramping.png'.format(tag), dpi=300, transparent=True)
    #plt.close(fig)
    
    #fig, ax = plt.subplots(1,1,figsize=(6,3))
    df_model_cs_ratio(df)  
    '''
    
    #plt.show()
    
    

def main_panel_fig_cs_intro():
    
    if True:
        id_mapper   =   use.create_id_mapper()
        for model in range(54):     #[2,5,10,15,27,28,33,39,47]:    # range(40):
            model = id_mapper[model]
            find_example_traces(model)
    if False:
        find_example_trace2_single()
    
    # morphology; set to true to plot
    if False:
        plot_morphology_with_ramping()
    
    # example traces; set to true to plot
    if False:
        f2,  a2 = plt.subplots(2,1, figsize=(5,6), gridspec_kw = {'height_ratios':[4,1], 'hspace':0})
        plot_example_traces(15, ax=a2[0])
        plot_alpha_func(ax=a2[1])
        for a in a2: 
            a.set_xlim([minT-50,maxT])
        a2[0].set_ylim([-90,40])
        f2.savefig('Figures/cs_example_ramping.png', dpi=300, transparent=True)
        #plt.show()
    
    # spike prob and time 2 spike over conditions
    if False:
        f1,a1       = plt.subplots(1,1, figsize=(6,3))
        f2,a2       = plt.subplots(1,1, figsize=(6,3))
        spike_array = extract_spikes(load=1)
        index       = 0 
        prevind     = 0
        nbg         = 5
        nmodels     = 54
        nbins       = len(conditions)
        COL         = ['k','#fdc086','#7fc97f','#beaed4']
        halfwidth   = 0.5
        
        spike_height = []
        # loop over conditions
        for i,ii in enumerate([1,2,10,20]):
            index = index + nbg*ii*nmodels
            
            # get only lines corresponding to one condition
            sub_array   = spike_array[prevind:index,:]
            
            # calculate (pspike and time to spike)
            first_column = sub_array[:,0]
            pspike = np.count_nonzero(first_column)/len(first_column)
            print(len(first_column), first_column.shape[0], pspike)
            t2s = np.mean( first_column[ first_column > 0 ] ) - 1000
            print('\t', t2s)
            
            spike_height.append(pspike)
            prevind = index
            
            # pspike
            a1.fill_between([i-halfwidth,i+halfwidth],[pspike,pspike], color=COL[i])
            # time to spike
            a2.fill_between([i-halfwidth,i+halfwidth],[t2s,t2s], [200,200], color=COL[i])
            
        
        a1.plot([-0.5,3.5],[0,0], 'k', lw=2)
        a1.plot([-0.5,3.5],[0.5,0.5], 'lightgrey')
        a1.plot([-0.5,3.5],[1,1], 'lightgrey')
        a1.axis('off')
        f1.savefig('Figures/ramping_pspike_hist.png', dpi=300, transparent=True)
        a2.plot([-0.5,3.5],[160,160], 'k', lw=2)
        a2.plot([-0.5,3.5],[200,200], 'lightgrey')
        a2.plot([-0.5,3.5],[300,300], 'lightgrey')
        a2.plot([-0.5,3.5],[400,400], 'lightgrey')
        a2.axis('off')
        f2.savefig('Figures/ramping_t2s_hist.png', dpi=300, transparent=True)
        plt.show()
            
    
        
    # bg (raster); set to true to plot
    if False:
        with open('Results/single_inVivo_ramping_{}_model27.json'.format('ACh+DA'), 'r') as handle:
            data = json.load(handle)
        
        fig = plt.figure(figsize=(5,3))
        grid    = plt.GridSpec(6, 6, hspace=0, wspace=0)
        ahist   = fig.add_subplot(grid[:1,:-1])
        ax      = fig.add_subplot(grid[1:,:-1])
        plot_bg_raster(data, ax=ax, ahist=ahist)
        ax.set_yticks([])
        ax.set_xlim([-505,505])
        ahist.set_xlim([-505,505])
        ahist.axis('off')
        fig.savefig('Figures/raster_background_example.png', dpi=300, transparent=True)
    
    

def cs_correlation_with_naf(load=False, subset='UniformMod'):
    
    # function also exists in extract_data.py
    
    
    #t = 'mixPat'
    t = 'ramping'
    
    ranges  = {'naf':(0.7,1.0), 'kaf':(0.7,1.0)}
    
        
    if load:
        with open('Results/Naf_factor_correlation/spike_probability_{}.json'.format(subset), 'r') as handle:
            data = json.load(handle)
        cs = data['complex']
        spikes = data['regular']
    else:
        cs = {  'naf':{'1':[], '0':[]},
                'kaf':{'1':[], '0':[]},
                'kir':{'1':[], '0':[]}
                }
        spikes = {  'naf':{'1':[], '0':[]},
                    'kaf':{'1':[], '0':[]},
                    'kir':{'1':[], '0':[]}
                    }
        for channel in ranges:
            files   = glob.glob('Results/Naf_factor_correlation/inVivo_ramping_vsFactor{}-{}*.json'.format(subset,channel))
            for File in files:
                try:
                    with open(File, 'r') as f: 
                        D = json.load(f)
                except:
                    print(f)
                    continue
                
                time = D['0']['time']
                
                for data in D.values():
                    
                    for d in data['data'].values(): 
                        # check trace for spikes (regular and complex)
                        spike_train = use.getSpikedata_x_y(time,d)
                        b, indx, trace = f4a.check_sliding_average_lowSampleData(d, 
                                                                    return_cs_index=1, 
                                                                    threshold=-37)
                        if b: 
                            w = 1
                            #plt.plot(time,d)
                            #plt.title('{}, {}'.format(channel, data['meta']['factors'][channel]))
                            #plt.show()
                        else: w = 0
                        
                        cs[channel][str(w)].append( data['meta']['fDA']['f'][channel] )
                        spikes[channel][str(w)].append(len(spike_train))
                    
                
        # save to file
        with open('Results/Naf_factor_correlation/spike_probability_{}.json'.format(subset), 'wt') as handle:
            json.dump({'regular':spikes,'complex':cs}, handle, indent=4)
    
    for chan in ranges:                            
        plt.figure(figsize=(4,2))
        plt.hist(cs[chan]['0']+cs[chan]['1'], bins=30, color='grey', histtype='stepfilled')
        plt.hist(cs[chan]['1'],               bins=30, color='r',    histtype='stepfilled', alpha=0.5)  
        #plt.title(chan)
        plt.axis('off')
        plt.savefig('Figures/cs_{}_dependence_ramping_{}.png'.format(chan,subset), dpi=300, transparent=True)
    
    #plt.figure()
    #plt.hist(spikes[0], color='k', range=(0,14), histtype='step', alpha=0.5)
    #plt.hist(spikes[1], color='r', range=(0,14), histtype='step', alpha=0.5)
    
    plt.show()
    
        
    
    #return cs

 

# -------------------------------------------------------------------------------------
# 2nd level

def plot_bg_raster(data, ax=None, ahist=None):
    
    if not ax:      fig,ax  = plt.subplots(1,1)
    if not ahist:   f,ahist = plt.subplots(1,1)
    
    # plot 
    color = {'glut':'k','gaba':'grey'}
    SPIKES = {'glut':[],'gaba':[]}
    for stype,a in zip(['glut','gaba'],[0,len(data['spikes']['ramp']['glut'])+10]):
        for j,key in enumerate(data['spikes']['ramp']['glut'].keys()):
            spikes = data['spikes']['ramp' ][stype][key] + \
                     data['spikes']['const'][stype][key]
            y = [j+a for s in spikes]
            ax.plot(spikes, y, '.', ms=3, c=color[stype])
            SPIKES[stype] += spikes
    
        ahist.hist(SPIKES[stype], range=(-500,500), bins=30, histtype='step', color=color[stype])
            
        
def plot_morphology_with_ramping(ax=None):
    
    if not ax: fig,ax = plt.subplots(1,1)
    
    import morph_lib_creator
    morph_with_none, sec_coordinates, stem2plot, sec2stem, morphology = morph_lib_creator.create(
                swc_file='../../morphologies/WT-dMSN_P270-20_1.02_SGA1-m24.swc')
    x = np.array(morph_with_none['x']).astype(float)
    y = np.array(morph_with_none['y']).astype(float)
    z = np.array(morph_with_none['z']).astype(float)
    
    # plot dendrites and soma
    ax.plot(x,y, color='grey')
    ax.plot(0,0,'o',color='grey', ms=10)
    
    # plot approximate synaptic locations
    for i,s in enumerate(sec_coordinates):
        # get coordinates
        coordinates = sec_coordinates[s]
        if i%2 == 0: c = 'k'
        else: c = 'grey'
        xi = coordinates[0]
        yi = coordinates[1]
        zi = coordinates[2] 
        for a,j in enumerate([24,16,8]):
            ax.plot(xi,yi, c=c, marker='.', ls='', ms=j, alpha=(a+1)/4, mew=0)
    ax.axis('off')
    fig.savefig('Figures/morphology_with_synaptic_bombardment.png', 
                dpi=300, 
                transparent=True)  
    plt.show()  
    

 

def inhibition_fig():
    
    # plot results from inhibition_ramping(_GABA).py
    
    cell_index = 10
    
    colors = {  'soma':'k',
                'dend[3]':'k',
                'dend[13]':'#525252',
                'dend[30]':'#fc4e2a',
                'dend[24]':'#e31a1c',
                'ctrl': 'k',
                'slow':'#8c510a',
                'fast':'#dfc27d'
                }
    maxT = 500;     minT =  100
    maxI = 0.0005;  minI = -0.035
    maxV = 40;      minV = -80
    alpha = 0.5;    ls   =  '--'
    compartments         =  ['soma', 'dend[24]', 'dend[30]'] #'dend[8]'
    
    conditions = [  {'type':'GABA-', 'dur':'', 'kin':'fast', 'comp':['soma', 'dend[24]', 'dend[30]']},
                    {'type':'GABA-', 'dur':'', 'kin':'slow', 'comp':['soma', 'dend[24]', 'dend[30]']},
                    {'type':'', 'dur':'', 'kin':''},
                    {'type':'', 'dur':'_short', 'kin':''}]
    
    # first panel (slow and fast gaba) 
    if False:                  
        for sim in conditions[:2]:
            it = sim['type']
            dur = sim['dur']
            kin = sim['kin']
            
            with open('Results/local_inh{}CS_ACh+DA_model{}{}{}.json'.format(it, cell_index, dur, kin), 'r') as f:
                data = json.load(f)
        
            time    = data['ctrl']['time']
            ctrl    = data['ctrl']
            
            for sec in sim['comp']:
                
                fig,ax = plt.subplots(2,1, figsize=(4,5), gridspec_kw = {'height_ratios':[5,2],'hspace':0})
                
                ax[0].plot([minT,maxT], [-63,-63],   '--k',       lw=2)
                ax[0].plot([minT,maxT], [-25,-25],   'lightgrey', lw=1)
                ax[0].plot([minT,maxT], [maxV,maxV], 'lightgrey', lw=1)
                
                # vm soma
                ax[0].plot(time, ctrl['vm']['soma'], color=colors['ctrl'], ls=ls, alpha=alpha)
                ax[0].plot(time, data[sec]['vm']['soma'], color=colors[kin])
                
                
                ax[0].set_ylim([minV+1,maxV+1])
                
                # inhibition
                ax[1].plot([0,500], [0,0], color=colors['ctrl'], ls=ls, alpha=alpha, lw=2)
                ax[1].plot(time, np.array(data[sec]['inh']), color=colors[kin], lw=2)
                ax[1].set_ylim([-0.01,0.1])
                if kin == '_fast':
                    ax[1].plot([400,450], [0.05,0.05], 'k', dash_capstyle='butt', lw=2)
                    ax[1].plot([400,400], [0.05,0.03], 'k', dash_capstyle='butt', lw=2)
                    
                for i,a in enumerate(ax):
                    a.set_xlim([minT,maxT])
                    a.axis('off')
                fig.savefig('Figures/local_inh{}CS_{}{}{}.png'.format(it,sec,dur,kin), dpi=300, transparent=True)
    
    
    # agregated results
    # results from inh_besk (candidates with cs + inhibition)
    def cosmetics(ax): 
        ax.plot([x[0],0], [0,0],        'k',         lw=2, ls='-', zorder=0)
        ax.plot([x[0],x[0]], [0,-0.02], 'k',         lw=2, ls='-', zorder=0)
        ax.plot([-30,-30], [0,-0.02],   'k',         lw=2, ls='-', zorder=0)
        ax.plot([x[0],0], [0.5,0.5],    'lightgrey', lw=1, ls='-', zorder=0)
        ax.plot([x[0],0], [1,1],        'lightgrey', lw=1, ls='-', zorder=0)
    def plot(ax,cond='fast'):
        if l == 'soma': 
            c = {'slow':colors['slow'], 'fast':colors['fast']} 
            ls = '-'
        else: 
            c = {'slow':colors[l], 'fast':colors[l]}   
            ls = '--' 
        ax.plot(x, both[cond], c[cond], marker='o', ls=ls)
    def ranges(ax):
        ax.set_ylim([-0.1,1.1])
        ax.axis('off')
    
    if True:
        
        dt          = [0,10,20,30,50,70,100,200]
        location    = ['3', '13', '24', '30', '']
        
        
        with open('../all_cs_candidates_ramping.json', 'r') as handle:
            candidates = json.load(handle)
        
        count = get_inhibit_data(load=True)
        
        f1,a1 = plt.subplots(1,1, figsize=(5,4))
        f2,a2 = plt.subplots(1,1, figsize=(3,3))
        f3,a3 = plt.subplots(1,1, figsize=(3,3))
        
        for loc in location:
            
            tag = '_2x' # '' for 1 slow and 10 fast (1x)
            
            if loc == '':                
                x=[]; slow=[]; fast=[]; serr = []; ferr = []
                for d in reversed(dt):
                    f=[]; s=[]
                    for bg in range(5):
                        f.append( sum(count[loc][tag][str(d)][str(bg)]['fast']) / sum(count[loc][tag][str(d)][str(bg)]['ctrl']) )
                        s.append( sum(count[loc][tag][str(d)][str(bg)]['slow']) / sum(count[loc][tag][str(d)][str(bg)]['ctrl']) )
                    fast.append( np.mean(f) )
                    slow.append( np.mean(s) )
                    ferr.append( np.std(f) )
                    serr.append( np.std(s) )
                    x.append(-1*d )
                both = {'fast':np.array(fast), 'slow':np.array(slow)}
                Mean = {'fast':np.array(fast), 'slow':np.array(slow)}
                Err  = {'fast':np.array(ferr), 'slow':np.array(serr)}
                l='soma'
                plot(a1,'fast')
                plot(a1,'slow')
            
            x=[]; slow=[]; fast=[]
            for d in reversed(dt):
                # only 1 bg (same as illustrated with NMDA; bg=1)
                fast.append( sum(count[loc][tag][str(d)][str(1)]['fast']) )
                slow.append( sum(count[loc][tag][str(d)][str(1)]['slow']) )
                x.append(-1*d )
            norm = sum(count[loc][tag]['0'][str(1)]['ctrl'])
            both = {'fast':np.array(fast)/norm, 'slow':np.array(slow)/norm}
            if loc == '':
                l='soma'
            else:
                l = 'dend[{}]'.format(loc)
            plot(a2,'fast')
            plot(a3,'slow')
        
        for cond in ['fast','slow']:
            a1.fill_between(x, Mean[cond]+Err[cond], Mean[cond]-Err[cond], color=colors[cond], alpha=0.1, zorder=0, lw=0)
        
        for a in [a1,a2,a3]: cosmetics(a)    
        for a in [a1,a2,a3]: ranges(a)
        
        f1.savefig('Figures/local_inh_CS_dt{}_soma.png'.format(tag), dpi=300, transparent=True)
        f2.savefig('Figures/local_inh_CS_dt{}_dend-fast.png'.format(tag), dpi=300, transparent=True)
        f3.savefig('Figures/local_inh_CS_dt{}_dend-slow.png'.format(tag), dpi=300, transparent=True)
        
        #plt.close('all')
        plt.show()
    
    # example nmda
    if False:
        
        with open('Results/local_inhCS_ACh+DA_model10.json', 'r') as f:
            data = json.load(f)
        
        D       = data['ctrl']
        time    = D['time']
        index   = next(i for i,t in enumerate(time) if t > 310)
        m1 = -0.02; m2 = -0.02; m3 = -0.02
        
        fig,ax = plt.subplots(1,1, figsize=(4,3))
        
        for comp in data['soma']['vm']:
            if comp == 'soma': continue
            # dendritic nmda and vm
            if comp in ['soma', 'dend[24]', 'dend[30]', 'dend[3]', 'dend[13]']: 
                #ax[0].plot(time, D['vm'][comp],   color=colors[comp], zorder=1e6)
                ax.plot(time, D['nmda'][comp], color=colors[comp], zorder=1e6)
            else: 
                #ax[0].plot(time, D['vm'][comp],   color='lightgrey', alpha=alpha)
                ax.plot(time, D['nmda'][comp], color='lightgrey', alpha=alpha)
            
            # min value
            '''
            if D['nmda'][comp][index] > m1:
                m1 = D['nmda'][comp][index]
                sec1 = comp
            elif D['nmda'][comp][index] > m2:
                m2 = D['nmda'][comp][index]
                sec2 = comp
            elif D['nmda'][comp][index] > m3:
                m3 = D['nmda'][comp][index]
                sec3 = comp'''
        
        # plot min
        '''    
        for sec in [sec1,sec2,sec3]:    
            ax[0].plot(time, D['vm'][sec],   color='k')
            ax[1].plot(time, D['nmda'][sec], color='k')
        print(sec1,sec2,sec3) '''
        
        ax.plot([minT,maxT],[0,0], '--k', lw=2) 
        ax.set_ylim([-0.04,0.001])
        ax.set_xlim([minT,maxT])
        ax.axis('off')
        fig.savefig('Figures/local_inhGABA-CS_nmda-example.png', dpi=300, transparent=True)
        plt.show()
        

  
def get_inhibit_data(load=False):
    
    if load:
        with open('Results/cs_following_inh.json', 'r') as handle:
            count = json.load(handle)
    else:
        dt          = [0,10,20,30,50,70,100,200]
        ninh        = ['_2x'] # '', '_2x'
        location    = ['3', '13', '24', '30', '']
        
        with open('../all_cs_candidates_ramping.json', 'r') as handle:
            candidates = json.load(handle)
            
        count   = {}
        for loc in location:
            count[loc] = {}
            for tag in ninh:
                fig,ax = plt.subplots(1,1, figsize=(5,4))
                print(tag)
                count[loc][tag] = {}
                for d in dt:
                    count[loc][tag][str(d)] = {}
                    for bg in range(5):
                        count[loc][tag][str(d)][str(bg)] = {    'fast':[],
                                                                'slow':[],
                                                                'ctrl':[],
                                                                'f':0   }
                    if not loc == '': 
                        l1 = '_dend[[]{}[]]'.format(loc)
                    else:
                        l1 = loc
                        
                    file_string = 'Results/Inhibition_of_CS/inVivo_ramping+inh_ACh+DA_candidate*_dt-{}{}{}.json'.format(d,tag,l1)
                    files = glob.glob( file_string )
                    print(file_string)
                    print(len(files))
                    
                    for f in files:
                        with open(f,'r') as handle:
                            data = json.load(handle)
                        
                        # get bgid
                        cand = f.split('candidate')[1].split('_')[0]
                        bg = candidates[int(cand)]['bg']
                        
                        if 'failed' in data['ctrl']: 
                            count[loc][tag][str(d)][str(bg)]['f'] += 1
                            continue
                        
                        for cond in ['ctrl', 'slow', 'fast']:
                            b = f4a.check_sliding_average_lowSampleData(data[cond]['vm'], threshold=-37)
                            if b: w = 1
                            else: w = 0
                            count[loc][tag][str(d)][str(bg)][cond].append(w)
            
        with open('Results/cs_following_inh.json', 'w') as handle:
            json.dump(count, handle, indent=4)
    
    return count                        
        
        
        
    
              
         
    

def plot_single_sodium(data, ax=None):
    
    if not ax: fig,ax = plt.subplots(1,1)
    
    #area = {'adist':94.24777960769379,'aprox':94.24777960769379,'soma':467.5946359395505}
    
    for key,trace in data['na']['I'].items():
        ax.plot(data['time'], trace, c=csingle[key])
    #for stype in ['soma','prox','dist']:
    #    trace = data['ACh+DA']['current'][4]['intr'][stype]['naf']
    #    ax.plot(data['ACh+DA']['time'], trace, c=c[stype])  
    
    ax.plot([minT,maxT], [0,0], '--k', lw=2, zorder=1e6)
    ax.plot([100,150], [-1.0,-1.0], '-k', lw=2)
    ax.plot([100,100], [-0.5,-1.0], '-k', lw=2)
    ax.plot([minT,maxT],[-2.2,-2.2],'lightgrey', zorder=1)
    #ax.plot([0,0],[-2.2,0],'k', lw=2, zorder=1)
    #ax.fill_between([0,16],[-2.0,-2.0],[-2.18,-2.18], color='grey', alpha=0.3, zorder=0,lw=0)
    ax.set_xlim(minT,maxT+1)
    ax.set_ylim(-2.5,0.1)
    
    ax.axis('off')


def plot_single_ca(data, ax=None):
    
    if not ax: fig,ax   = plt.subplot(1,1)
    
    tm = data['time']
    
    colCa   = {'ca':'#fdb863','cal':'#5e3c99'}
    
    cal     = data['Ca']['LT-NMDA']
    ca      = data['Ca']['NPQR']
    acal    = np.zeros((len(cal),len(tm)))
    aca     = np.zeros((len(ca),len(tm)))
    max_cal = 0
    for i,key in enumerate(ca.keys()):
        trace = ca[key]
        aca[i,:] = trace
        ax.plot(tm,trace, c=colCa['ca'], alpha=0.4, zorder=i+1e3)
        trace = cal[key]
        acal[i,:] = trace
        ax.plot(tm,trace, c=colCa['cal'], alpha=0.2)
        if max(trace) > max_cal: max_cal = max(trace); dkey = key 
    mcal = np.mean(acal, axis=0)
    mca  = np.mean(aca , axis=0)
    ax.plot(tm,mcal, c=colCa['cal'], alpha=1.0, lw=2, zorder=i+1+1e3)
    ax.plot(tm,mca,  c=colCa['ca'], alpha=1.0, lw=2, zorder=i+1+1e3)
    x1,x2,y1,y2=ax.axis()
    ax.set_ylim([-0.01,0.201])
    ax.plot([minT,maxT],[0.2,0.2],'lightgrey', zorder=1)
    ax.plot([minT,maxT],[0,0],'--k',lw=2)
    ax.set_xlim(minT,maxT+1)
    ax.axis('off')
    
    #print(dkey, data['dist'][dkey], [data['Vm']['dendrites'][dkey][i] for i,t in enumerate(data['time']) if t > 330][0])



def plot_single_glut(data, ax=None): 
    
    if not ax: fig,ax   = plt.subplot(1,1)
    
    tm = data['time']
    
    # ampa
    ampa = data['ampa']
    array = np.zeros((len(ampa),len(tm)))
    for i,key in enumerate(ampa.keys()):
        trace = ampa[key]
        ramp = data['ramp']['ampa'][key]
        combined = np.add(trace,ramp)
        array[i,:] = combined
        ax.plot(tm,trace, 'lightgreen')
    mampa = np.mean(array, axis=0)
    
    # nmda
    nmda = data['nmda']['I']
    max_nmda = 0
    for i,key in enumerate(nmda.keys()):
        trace = nmda[key]
        ramp = data['ramp']['I'][key]
        combined = np.add(trace,ramp)
        array[i,:] = combined
        ax.plot(tm,trace, 'lightgrey')
        #if min(trace) < max_nmda: max_nmda = min(trace); dkey = key 
    marray = np.mean(array, axis=0)
    #ind = [i for i,t in enumerate(tm) if t > 310] # what is these fore?
    #argmin = np.argmin(marray[ind])    # what is these fore?
    ax.plot(tm,mampa, 'green')
    ax.plot(tm,marray, 'grey')
    ax.set_ylim([-0.05,0.001])
    ax.set_xlim(minT,maxT+1)
    ax.plot([minT,maxT],[0,0],'--k',lw=2)
    ax.axis('off')
    
    ax.plot([minT,maxT],[-0.045,-0.045],'lightgrey', zorder=1)
    #ax.plot([0,0],[-0.045,0],'k', lw=2)
    #ax.fill_between([0,16],[-0.045,-0.045],[-0.048,-0.048], color='grey', alpha=0.3, zorder=0,lw=0)
    
    #ax.plot([150,200], [-0.043,-0.043], '-k', lw=2)
    ax.plot([minT+1,minT+1], [-0.033,-0.043], '-k', lw=2)    


def correlate_nmdaDend_with_vmSoma(data, ax=[]):
    # TODO: clean up
    
    if not len(ax): fig,ax   = plt.subplots(2,1, figsize=(4,4))  
    from scipy import stats   
    
    tm = data['time']
    
    # correlation between mean nmda current and 
    #                     low pass filtered vm
    array = np.zeros((len(data['nmda']['I']),len(tm)))
    aampa = np.zeros((len(data['ampa']),len(tm)))
    
    for i,key in enumerate(data['nmda']['I'].keys()):
        trace = data['nmda']['I'][key]
        ramp  = data['ramp']['I'][key]
        combined = np.add(trace,ramp)
        array[i,:] = combined
        trace = data['ampa'][key]
        ramp = data[ 'ramp']['ampa'][key]
        combined = np.add(trace,ramp)
        aampa[i,:] = combined
    
    nmda = np.mean( array, axis=0 )
    ampa = np.mean( aampa, axis=0 )
        
    b, index, smooth = f4a.check_sliding_average_lowSampleData(data['Vm']['soma'], return_cs_index=True, threshold=-40)
    
    ax[0].plot([0,0], [-1,1], '--k', lw=2)
    # normalize traces (to fit in same fig)
    p = {   0:{'c':'lightgreen', 'shift':-1, 'a':1},
            1:{'c':'grey', 'shift':-1, 'a':1},
            2:{'c':'k', 'shift':0, 'a':0.1}}
    for i,trace in enumerate([ampa,nmda,data['Vm']['soma']]):
        mi = min(trace)
        ma = max(trace)
        normed =  (np.array(trace)-mi)/(ma-mi)+p[i]['shift']
        ax[0].plot(tm, normed, c=p[i]['c'], alpha=p[i]['a'])
    normed = (np.array(smooth)-mi)/(ma-mi)
    ax[0].plot(tm, normed, c='k')
    
    # 2:{'ls':'-','c':'k'}
    corr = []; campa = []
    corr.append( stats.pearsonr(smooth, nmda)[0] )
    campa.append( stats.pearsonr(smooth, ampa)[0] )
    print('***** CORRELATION ALL TRACE **********')
    print( 'NMDA',corr[0],'\tAMPA', campa[0] )
    print('------ CORRELATION FIRST SEC ---------')
    T = 0
    PARTIAL_SMOOTH = [smooth[i] for i,t in enumerate(tm) if t < T]
    PARTIAL_NMDA = [nmda[i] for i,t in enumerate(tm) if t < T]
    PARTIAL_AMPA = [ampa[i] for i,t in enumerate(tm) if t < T]
    corr.append( stats.pearsonr(PARTIAL_SMOOTH, PARTIAL_NMDA)[0] )
    campa.append( stats.pearsonr(PARTIAL_SMOOTH, PARTIAL_AMPA)[0] )
    print( 'NMDA',corr[1],'\tAMPA', campa[1] )
    print('------ CORRELATION RAMP ---------')
    PARTIAL_SMOOTH = [smooth[i] for i,t in enumerate(tm) if t > T]
    PARTIAL_NMDA = [nmda[i] for i,t in enumerate(tm) if t > T]
    PARTIAL_AMPA = [ampa[i] for i,t in enumerate(tm) if t > T]
    corr.append( stats.pearsonr(PARTIAL_SMOOTH, PARTIAL_NMDA)[0] )
    campa.append( stats.pearsonr(PARTIAL_SMOOTH, PARTIAL_AMPA)[0] )
    print( 'NMDA',corr[2],'\tAMPA', campa[2] )
    print('***** END CORRELATION **********')
    
    width = 0.35
    x = np.arange(len(corr))
    ax[1].bar(x - width/2, corr,  width, color=p[1]['c'])
    ax[1].bar(x + width/2, campa, width, color=p[0]['c'])
    
    ax[1].plot([-0.5,2.5], [0,0], '--k', lw=2)
    ax[1].plot([-0.5,2.5], [-1,-1], 'lightgrey')
    ax[1].set_ylim([-1.1, 0.1])
    for a in ax:
        a.axis('off')


def plot_single_vm_diff(data, ax=[]):
    
    if not len(ax): fig,ax = plt.subplot(1,1)
    
    # plot
    tm      = data['time']
    adist   = data['Vm']['axon']['adist']
    aprox   = data['Vm']['axon']['aprox']
    soma    = data['Vm']['soma']  
    
    ax[0].fill_between(tm, adist, aprox, color=csingle['adist'], lw=0)
    ax[0].fill_between(tm, aprox, soma,  color=csingle['aprox'], lw=0)
    ax[0].plot([minT,maxT],[-25,-25],'lightgrey')
    ax[0].plot([minT,maxT], [-60,-60], '--k', lw=2, zorder=1)
    for sec in [soma, aprox, adist]:
        ax[0].plot(tm, sec, 'k')
    ax[0].set_ylim([-61,-15])
    
    ax[1].plot([minT,maxT],[8,8],'lightgrey')
    ax[1].plot([minT,maxT], [0,0], '--k', lw=2, zorder=1)
    ax[1].plot(tm, np.array(adist)-aprox, 'r')
    ax[1].plot(tm, np.array(aprox)-soma,  'm')
    ax[1].set_ylim(-1,8.1)
    for a in ax:
        a.set_xlim([310,400])
        a.axis('off')

def plot_single_vm_sodium_zoomed(data, ax=[]):
    
    if not len(ax): fig,ax = plt.subplot(1,1)
    
    # plot
    tm      = data['time']
    adist   = data['Vm']['axon']['adist']
    aprox   = data['Vm']['axon']['aprox']
    soma    = data['Vm']['soma']  
    
    #ax[0].fill_between(tm, adist, aprox, color=csingle['adist'], lw=0)
    #ax[0].fill_between(tm, aprox, soma,  color=csingle['aprox'], lw=0)
    ax[0].plot([minT,maxT],[-25,-25],'lightgrey')
    ax[0].plot([minT,maxT], [-60,-60], '--k', lw=2, zorder=1)
    for sec,key in zip([soma, aprox, adist],['soma', 'aprox', 'adist']):
        ax[0].plot(tm, sec, c=csingle[key])
    ax[0].set_ylim([-61,-15])
    
    #ax[1].plot([minT,maxT],[8,8],'lightgrey')
    ax[1].plot([minT,maxT], [0,0], '--k', lw=2, zorder=1)
    for key,trace in data['na']['I'].items():
        ax[1].plot(tm, trace, c=csingle[key])
    ax[1].set_ylim(-0.3,0.01)
    for a in ax:
        a.set_xlim([310,400])
        a.axis('off')
    

def plot_single_vm_all(data, ax=[]): 
    # remove?
    if not len(ax):  fig,ax  = plt.subplot(1,1)
    
    ax[0].plot([minT,maxT], [-63,-63], '--k', lw=2, zorder=1e6)
    ax[0].plot([maxT,maxT], [-63, 40], 'lightgrey', lw=1.5, zorder=1e6)
    #ax[0].fill_between([0,16],[-63,-63],[-70,-70], color='grey', alpha=0.3, zorder=0,lw=0)
    
    lw= 2
        
    tm = data['time']
    array = np.zeros((len(data['Vm']['dendrites']),len(tm)))
    for i,trace in enumerate(data['Vm']['dendrites'].values()):
        array[i,:] = trace
        ax[0].plot(tm,trace, 'lightgrey', alpha=0.3)
    marray = np.mean(array, axis=0)
    ax[0].plot(tm,data['Vm']['axon']['aprox'], c=csingle['aprox']) 
    ax[0].plot(tm,data['Vm']['axon']['adist'], c=csingle['adist'])   
    ax[0].plot(tm,data['Vm']['soma'], c=csingle['soma'], lw=lw)
    ax[0].plot(tm,marray,  'grey', lw=lw)
    
    maxV = f4a.plot_window_current(ax[1], color='r', minV=-63, return_peak=True)
    ax[0].plot([minT,maxT],[maxV,maxV],'lightgrey', zorder=0)
    ax[0].plot([minT,maxT],[40,40],'lightgrey')
    
    ax[0].axis('off')
    ax[1].axis('off')
    
    ax[0].set_ylim(-80,41)
    ax[0].set_xlim(minT,maxT+1)
    ax[1].set_ylim(-80,41)
    ax[1].set_xlim(-0.01,1.02)
    
    print(maxV)
    
    return marray

def plot_single_vm_soma_axon(data, ax=[]): 
    
    if not len(ax):  fig,ax  = plt.subplot(1,1)
    
    ax[0].plot([minT,maxT], [-63,-63], '--k', lw=2, zorder=1e6)
    ax[0].plot([maxT,maxT], [-63, 40], 'lightgrey', lw=1.5, zorder=1e6)
    
    lw= 2
        
    tm = data['time']
    ax[0].plot(tm,data['Vm']['axon']['aprox'], c=csingle['aprox']) 
    ax[0].plot(tm,data['Vm']['axon']['adist'], c=csingle['adist'])  
    ax[0].plot(tm,data['Vm']['soma'], c=csingle['soma'], lw=lw)
    # na window
    maxV = f4a.plot_window_current(ax[1], color='r', minV=-63, return_peak=True)
    ax[0].plot([minT,maxT],[maxV,maxV],'lightgrey', zorder=0)
    ax[0].plot([minT,maxT],[40,40],'lightgrey')
    #ax[0].plot([0,0],[-63,40],'k', lw=2)
    ax[0].set_xlim(minT,maxT+1)
    ax[1].set_xlim(-0.01,1.02)
    for a in ax:
        a.axis('off')
        a.set_ylim(-80,41)
    
    print(maxV)



def plot_single_vm_soma_dend(data, ax=[]): 
    
    if not len(ax):  fig,ax  = plt.subplot(1,1)
    
    ax[0].plot([minT,maxT], [-63,-63], '--k', lw=2, zorder=1e6)
    ax[0].plot([maxT,maxT], [-63, 40], 'lightgrey', lw=1.5, zorder=1e6)
    
    lw= 2
        
    tm = data['time']
    array = np.zeros((len(data['Vm']['dendrites']),len(tm)))
    for i,trace in enumerate(data['Vm']['dendrites'].values()):
        array[i,:] = trace
        ax[0].plot(tm,trace, 'lightgrey', alpha=0.3)
    marray = np.mean(array, axis=0)
    ax[0].plot(tm,data['Vm']['soma'], c=csingle['soma'], lw=lw)
    ax[0].plot(tm,marray,  c=csingle['dend'], lw=lw)
    # mg-block
    maxV = f4a.plot_mgblock(ax[1], c_mg='k', c_2nd='grey', p_1stDer=1, p_fbtw=0, return_peak=True, minV=-63)
    ax[0].plot([minT,maxT],[maxV,maxV],'lightgrey', zorder=0)
    ax[0].plot([minT,maxT],[40,40],'lightgrey')
    #ax[0].plot([0,0],[-63,40],'k', lw=2)
    
    ax[0].axis('off')
    ax[1].axis('off')
    
    ax[0].set_ylim(-80,41)
    ax[0].set_xlim(minT,maxT+1)
    ax[1].set_ylim(-80,41)
    ax[1].set_xlim(-0.01,1.02)
    
    print(maxV)
    
    return marray

def find_example_traces_ispn():
    
    #ax[0].set_title(mid)
    alpha = 0.5
    
    count = {2:0, 1:0, 0:0}
    
    flag = 0
    
    #files = glob.glob('Results/*_ACh+DA_*model[!0]*.json')
    #files = glob.glob('Results/*_redNaf*.json')
    #files = glob.glob('Results/*_incNMDA*.json')
    files = glob.glob('Results/*_ctrl*.json')
    fig,ax = plt.subplots(1,3, figsize=(12,4))
    for f in files:
        #f = 'Results/inVivo_ramping_{}_model{}.json'.format(c,mid)
        try:
            with open(f, 'r') as handle:
                data = json.load(handle)
            print(f)
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
                trace = data[i][str(bg)] #[data[str(i)]['0'][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                
                if max(trace) < -10:
                    ax[0].plot(time, trace, alpha=alpha)
                    count[0] += 1
                elif f4a.check_sliding_average_lowSampleData(trace, threshold=-37):
                    # check for complex spikes.
                    ax[2].plot(time, trace, alpha=alpha)
                    plt.plot(time, trace, alpha=alpha)
                    plt.show()
                    s = "'mod':{}, 'bg':{}".format(i,bg)
                    print('{'+s+'},')
                    count[2] += 1
                    flag=1
                else:
                    ax[1].plot(time, trace, alpha=alpha)
                    count[1] += 1
    
    print(count)            
    for i in range(3):
        ax[i].set_xlabel(count[i])
    if True:
        plt.show()
    else:
        plt.close('all')

def find_example_traces(mid, ax=None):
    
    if not ax: fig,ax = plt.subplots(1,3, figsize=(12,4))
    
    ax[0].set_title(mid)
    alpha = 0.5
    
    count = {2:0, 1:0, 0:0}
    
    flag = 0
    
    for c in ['ACh+DA']: # 'ctrl','ACh','DA',
        f = 'Results/inVivo_ramping_{}_model{}.json'.format(c,mid)
    
        with open(f, 'r') as handle:
            data = json.load(handle)
        
        tremove = 1000
        time = [t-1000 for t in data['time'] if t > tremove]
        # da fid: ach fid: condition: bg id
        for i in [0,2]:
            for j in range(10):
                for bg in data[str(i)][str(j)][c]:
                    trace = [data[str(i)][str(j)][c][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                    # check for APs and complex spikes.
                    b = f4a.check_sliding_average_lowSampleData(trace, threshold=-37)
                    if b:
                        ax[2].plot(time, trace, c=colors[c], alpha=alpha)
                        s = "'mod':{}, 'fach':{}, 'fda':{}, 'bg':{}".format(mid,i,j,bg)
                        print('{'+s+'},')
                        count[2] += 1
                        flag=1
                    elif max(trace) > -10:
                        ax[1].plot(time, trace, c=colors[c], alpha=alpha)
                        count[1] += 1
                    else:
                        ax[0].plot(time, trace, c=colors[c], alpha=alpha)
                        count[0] += 1
                    
    for i in range(3):
        ax[i].set_xlabel(count[i])
    if flag:
        plt.show()
    else:
        plt.close('all')

def find_example_trace2_single():
    
    candidates  = get_candidates()
    
    alpha       = 1.0
    tremove     = 1000
    
    #for i,candidate in enumerate(candidates):
    # [2,6,7,9,10,15,30,33,46,51,53,79,90,116]
    
    for i in [2,18,68]: #[2,6,7,9,10,15,30,33,46,51,53,79,90,116]: #[2,6,30,116]:  #[2,15,79]:
        
        candidate = candidates[i]
        
        fig,ax = plt.subplots(1,1)
        
        mid     = candidate['mod']
        fach    = candidate['fach']
        fda     = candidate['fda']
        bg      = candidate['bg']
        
        for c in ['ctrl','ACh','DA','ACh+DA']:
            f = 'Results/inVivo_ramping_{}_model{}.json'.format(c,mid)
        
            with open(f, 'r') as handle:
                data = json.load(handle)
            
            if 'DA'  in c:  
                da   = str(fda)
            else: da = '0'
            if 'ACh' in c: 
                ach   = str(fach)
            else: ach = '0'
            
            time = [t-1000 for t in data['time'] if t > tremove]
            trace = [data[ach][da][c][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
            
            ax.plot(time, trace, c=colors[c], alpha=alpha)
                    
        ax.set_title('cadidate {};   model:{}, ach:{}, da:{}, bg:{}'.format(i, mid,fach,fda,bg))
    plt.show()

def plot_cs_dipps(load=False):
    # TODO: move somewhere close to functions for detailed single trace (glut)
    
    fig, ax = plt.subplots(1,1,figsize=(4,4))
    
    dipps = get_cs_features(load=load)['dipps']
    
    ax.hist(dipps)
    ax.plot([-63,40], [0,0], '--k')
    ax.plot([-42,-42],[0,10], '--k')
    v = np.arange(-63,40)
    alpha=-0.062; beta=3.57
    mgblock =  1.0 / (1 + 1.0 * np.exp(alpha * v) / beta )
    diff = np.diff(mgblock)
    ac = np.diff(diff)
    ax.plot(v[1:-1],np.divide(ac,max(ac))*10, color='grey')
    ax.axis('off')
    fig.savefig('Figures/cs_initiation_zone_ramping.png', dpi=300, transparent=True)
    

def plot_cs_dur(load=False):
    # TODO: move somewhere close to functions for detailed single trace (glut)
    
    fig, ax = plt.subplots(1,1,figsize=(4,4))
    
    data = get_cs_features(load=load)
    dur = data['dur']
    mixtdur = get_mixed_features()
    
    n1,b1,p1 = ax.hist(dur, bins=20, range=(0,250))
    for p in p1:
        p.set_height(p.get_height()/len(dur)*200)
    n2,b2,p2 = ax.hist(mixtdur, bins=20, range=(0,250), color='lightblue', alpha=0.8)
    for p in p2:
        p.set_height(p.get_height()/len(mixtdur)*200)
    
    ax.plot([0,250], [0,0], 'k', lw=2)
    for x in [0,100,200]:
        ax.plot([x,x], [-2,0], 'k', lw=2)
    for x in [0,100,200]:
        ax.plot([x+50,x+50], [-1,0], 'k', lw=2)
    ax.axis('off')
    ax.set_ylim([-2,100])
    fig.savefig('Figures/cs_duration_ramping.png', dpi=300, transparent=True)
    plt.show()


def plot_multiple_nmdas_spikes_cs_dur():
    
    files = glob.glob('Results/single_inVivo_ramping_ACh+DA_model*_candidate*')
    for f in files:
        print(f)
        cid = f.split('candidate')[1].split('.')[0]
        with open(f, 'r') as handle:
            data = json.load(handle)  
        f,a = plt.subplots(2,1, figsize=(3,5), gridspec_kw = {'height_ratios':[1,2], 'hspace':0})  
        a[0].plot(data['time'],data['Vm']['soma'],'#1f77b4')
        for i,trace in enumerate(data['nmda']['I'].values()):
            if min(trace) < -0.018: c='#1f77b4'; alpha=1; o=1e6+i
            else: c='grey'; alpha=0.5; o=i
            a[1].plot(data['time'],trace, c=c, alpha=alpha, zorder=o)
        a[0].plot([330,480], [-60,-60], 'k')
        for ax in a:
            ax.set_xlim([100,500])
            ax.axis('off')
        a[1].set_ylim([-0.025,0])
        f.savefig('Figures/multiple_nmda-spikes_prolongs_CS_example_modl{}.png'.format(cid), dpi=300, transparent=False)



def get_mixed_features():
    with open('Results/extracted_cs_dipps_mixPat.json', 'r') as handle:
        features = json.load(handle) 
    
    return features['dur']    


def get_cs_features(load=False):
    
    if load:
        with open('Results/extracted_cs_dipps-1.json', 'r') as handle:
            features = json.load(handle)   
    else:
        dipps       = [] 
        durations   = []
        candidates  = get_candidates()
        tremove     = 1000
        for i,candidate in enumerate(candidates):    
            mid     = candidate['mod']
            fach    = candidate['fach']
            fda     = candidate['fda']
            bg      = candidate['bg']
            
            for c in ['ACh+DA']:
                f = 'Results/inVivo_ramping_{}_model{}.json'.format(c,mid)
            
                with open(f, 'r') as handle:
                    data = json.load(handle)
                
                if 'DA'  in c:  
                    da   = str(fda)
                else: da = '0'
                if 'ACh' in c: 
                    ach   = str(fach)
                else: ach = '0'
                
                time = [t-1000 for t in data['time'] if t > tremove]
                trace = [data[ach][da][c][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
                
                # get extreme values (peaks and dipps)
                b, indx, t = f4a.check_sliding_average_lowSampleData_from0(trace, 
                                                                            return_cs_index=1, 
                                                                            threshold=-37)
                try:
                    pads = f4a.extract_cs_dipps(trace,time, csi=indx[1])
                except:
                    print('FAILED! at dipp extraction')
                    continue
                #print(pads)
                mi = pads['min'][-1:]   # min index
                dipps.append(trace[mi[-1]])
                try:
                    dur = f4a.extract_cs_dur(trace,time, start_index=mi[-1])
                except:
                    print('FAILED!  at dur extraction')
                    continue
                print(dur)
                durations.append(dur)
                
                        
            #ax.set_title('cadidate {};   model:{}, ach:{}, da:{}, bg:{}'.format(i, mid,fach,fda,bg))
        
        # save dipps to file
        features = {'dipps':dipps, 'dur':durations}
        with open('Results/extracted_cs_dipps-1.json', 'w') as handle:
            json.dump(features, handle, indent=4)
        
    
    return features
    
    

def example_traces_factors():
    
    candidates  = get_candidates()
    
    alpha       = 1.0
    tremove     = 1000
    
    cid=10; rid=7
    
    candidate = candidates[cid]
    
    mid     = candidate['mod']
    fach    = candidate['fach']
    fda     = candidate['fda']
    bg      = candidate['bg']
    
    col     = {mid:'r', rid:'grey'}
    
    base    = 10
    dt      = 100
    
    for i,tag in zip([mid,rid],['cs','spikes']): 
        fig,ax = plt.subplots(1,1,figsize=(2,4)) 
        ax.plot([0,500], [-70,-70], c='k', lw=2, ls='--')
        f = 'Results/inVivo_ramping_ACh+DA_model{}.json'.format(i)
        with open(f, 'r') as handle:
            data = json.load(handle)
        da  = str(fda)
        ach = str(fach)
        
        time = [t-1000 for t in data['time'] if t > tremove]
        trace = [data[ach][da]['ACh+DA'][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
        
        ax.plot(time, trace, c=col[i], alpha=alpha)
        ax.set_xlim([0,500])
        ax.set_ylim([-80,40])
        if tag == 'cs':
            ax.plot([base,base+dt],[-20,-20],'k', lw=2)
            ax.plot([base,base],[0,-20],  'k', lw=2)
        plt.axis('off')
        fig.savefig('Figures/example_traces_{}_4factors_ramping.png'.format(tag), dpi=300, transparent=True)
    plt.show()

def plot_example_traces(cid, ax=None):
    
    if not ax: fig,ax = plt.subplots(1,1)
    
    candidates  = get_candidates()
    candidate   = candidates[cid]
    mid         = candidate['mod']
    fach        = candidate['fach']
    fda         = candidate['fda']
    bg          = candidate['bg']
    
    tremove     = 950
        
    for c in ['ctrl','DA','ACh','ACh+DA']:
        f = 'Results/inVivo_ramping_{}_model{}.json'.format(c,mid)
        
        with open(f, 'r') as handle:
            data = json.load(handle)
        
        if 'DA'  in c:  
            da   = str(fda)
        else: da = '0'
        if 'ACh' in c: 
            ach   = str(fach)
        else: ach = '0'
        
        time = [t-1000 for t in data['time'] if t > tremove]
        # da fid: ach fid: condition: bg id
        trace = [data[ach][da][c][str(bg)][ind] for ind,t in enumerate(data['time']) if t > tremove]
        
        if c == 'ACh+DA': lw = 2
        else:lw = 1.5
        ax.plot(time, trace, c=colors[c], lw=lw)
    
    lw=3
    ax.plot([minT-50,maxT], [-70, -70], '--k', lw=2, zorder=0) 
    
    xbase=150; ybase=-20; xlen=-100; ylen=-20
    ax.plot([xbase+xlen,xbase], [ybase,ybase], 'k', lw=lw)
    ax.plot([xbase,xbase], [ybase+ylen,ybase], 'k', lw=lw)
    ax.axis('off')
    #plt.show()


def plot_alpha_func(ax=None, tau=300, gmax=1):
    
    if not ax: fig,ax = plt.subplots(1,1)
    
    ax.plot([-50,500], [0.06,0.06], 'k', lw=2, ls='--')
    
    modOnTime   = 0
    tstop       = 500
    t           = np.arange(-50,tstop)
    
    transient = np.array([use.alpha(ht, modOnTime, 1, tau) if ht >= modOnTime else 0 for ht in t])
    
    for i,c in enumerate(['ACh','DA','ACh+DA']):
        ax.plot(t,transient+i*0.06, c=colors[c], lw=2)
    ax.axis('off')



# ==============================================================================


def df_proportion_condition(df):
    
    # filter frame to only use factors da in (0,2)
    df_filtered = df.loc[(df['fDA'] == 0) | (df['fDA'] == 2)]
    
    # get the number of models that contain complex spikes
    dcs = df_filtered.loc[(df_filtered['cs'] == 1)]
    print('::: number of models with cs = {}'.format( dcs.groupby('model')['cs'].nunique().count()))
    
    COL = ['grey', 'r']
    for i,cond in enumerate(['#spikes', 'cs']):
        fig, ax = plt.subplots(1,1, figsize=(4,2))
        df_pos      = df_filtered.loc[df_filtered[cond] >= 1]
        n1,b1,p1 = ax.hist(    df_filtered['condition'], 
                    range=(-0.5,3.5), 
                    bins=4, 
                    color='k', 
                    histtype='step',
                    lw=1, 
                    alpha=1.0
                    )
        n2,b2,p2 = ax.hist(    df_pos['condition'], 
                    range=(-0.5,3.5), 
                    bins=4, 
                    color=COL[i], 
                    histtype='stepfilled', 
                    alpha=1.0
                    )
        
        percent = n2/n1*100
        
        if False:
            # write percentage
            for x in range(0,4):
                ax.text(x+0.5*(i-1), n2[x]+100, '{0:.0f}%'.format(percent[x]), c=COL[i])
        
        print(cond, '-------------------------------')
        print('{}, {}'.format(n2,n1))
        print('%:', percent)
        
        ax.set_xticks([0,1,2,3])
        ax.set_xticklabels(conditions)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        fig.savefig('Figures/proportion_{}_condition_ramping.png'.format(cond.replace('#','')), dpi=300, transparent=True)
        
    
    #plt.show()


def df_proportion_model(df, ax=None, tag=''):
    
    if not ax: 
        fig,ax = plt.subplots(1,1, figsize=(6,3))
        f2, a2 = plt.subplots(1,1, figsize=(6,3))
    
    # ADDED FILTER (ONLY CLUST) HERE!!!
    df_pos      = df.loc[df['cs'] == 1]
    df_spikes   = df.loc[df['#spikes'] > 0]
    
    if tag == '':
        condition = 'ACh+DA'
    else:
        condition = 3
    
    # total number of each type (calculated from df_norm into normList)
    # 40 pattern x 3 c-lev x 5 bg = 600; 600 x Nmod (1,2,2,4)
    #condMap = {'ctrl':600,'ACh':1200,'DA':1200,'ACh+DA':2400}
    
    df_pos_cond    = df_pos.loc[   (df_pos[   'condition'] == condition) & ((df_pos['fDA']==0)|(df_pos['fDA']==2))]
    df_spikes_cond = df_spikes.loc[(df_spikes['condition'] == condition) & ((df_spikes['fDA']==0)|(df_spikes['fDA']==2))]
    
    # get number of trial of each model
    df_norm     = df.loc[(df['condition'] == condition) & ((df['fDA']==0)|(df['fDA']==2))]
    models  = np.array(sorted(df.model.unique()))
    cs_proportion   = []
    normList        = []
    spikes_prop     = []
    for m in models:
        cs_proportion.append(   len(df_pos_cond.loc[    df_pos_cond[    'model']==m]) )
        normList.append(        len(df_norm.loc[        df_norm[        'model']==m]) )
        spikes_prop.append(     len(df_spikes_cond.loc[ df_spikes_cond[ 'model']==m]) )
    
    # sort (on cs proportion--the spike plot is also identically sorted)
    sort_i_count    = np.argsort(cs_proportion)
    nbins           = len(models)
    theList         = []
    theSpikeList    = []
    for i,indx in enumerate(sort_i_count):
        theList += [i]*cs_proportion[indx]
        theSpikeList += [i]*spikes_prop[indx]
        
    # plot
    AX = [ax, a2]
    mmin = -0.5; mmax=54
    COL = ['r','grey']
    cs_proportion = []
    for i,a in enumerate(AX):
        
        if i == 0:
            L = theList
        else:
            L = theSpikeList
        n1,b1,p1 = a.hist( L, 
                    range=(-0.5,nbins-0.5),  
                    bins=nbins, 
                    color=COL[i],
                    histtype='bar', 
                    alpha=1.0,
                    ec='black'
                    )
        
        # normalize the bins (into percent cs)
        for j,item in enumerate(p1):
            height = item.get_height()
            item.set_height(        height/normList[j]*100 )
            if i == 0:
                cs_proportion.append(   height/normList[j]*100 )
        
        # sortera x-ticks (modelles) least -> most cs 
        a.set_xticks(range(len(n1)))
        a.set_xticklabels( models[sort_i_count], fontsize='small' )
        
        # finalize and save the plot window
        a.set_ylim([0,110])
        a.plot([mmin,mmax],[50,50],'lightgrey')
        a.plot([mmin,mmax],[100,100],'lightgrey')
        a.plot([mmin,mmax],[0,0],'k', lw=2)
        a.plot([mmax,mmax],[0,100],'k', lw=2)
        a.axis('off')
        a.invert_xaxis()
            
    print(models[sort_i_count])
    fig.savefig('Figures/cs_proportion-model-cond{}_ramping{}.png'.format(condition,tag), dpi=300, transparent=True)
    f2.savefig('Figures/nspikes_proportion-model-cond{}_ramping{}.png'.format(condition,tag), dpi=300, transparent=True)
    plt.show()
    
    # TODO: remove below if not used (EMBARGO). think rheobase is used at least.
    # correlate with model parameters ========================================================================================
    path            = '../Cholinergic_modulation'
    with open('{}/D1_71bestFit_updRheob.pkl'.format(path), 'rb') as f:
        model_sets  = pickle.load(f, encoding="latin1")
    
    # rheobase ----------------------------------------------------
    rheobase = []
    fig, ax = plt.subplots(1,1, figsize=(3,3)) #, figsize=(18,9) )
    
    # get rheobase for each model (sorted least to max cs)
    for indx in models[sort_i_count]:
        rheobase.append(model_sets[indx]['rheobase'])
    
    # plot rheobase vs cs proportion
    cs_pos_models = []
    for i in range(len(cs_proportion)):
        if cs_proportion[i] > 0:
            cs_pos_models.append(i)
            ax.plot(cs_proportion[i], rheobase[i], 'o', ms=10, mew=1, mec='w', c='k')
    plot_linreg(ax, np.array(cs_proportion)[cs_pos_models], np.array(rheobase)[cs_pos_models])
    fig.savefig('Figures/lincorr_csProp-rheobase-cond{}_ramping{}.png'.format(condition,tag), dpi=300, transparent=True)
    
    # model variables ----------------------------------------------------
    excDict     = {'naf':{}, 'kaf':{}}#, 'kas':{}, 'kir':{}, 'sk':{}} #, 'can':{}, 'c32':{}, 'c33':{}}
    listOfList  = []
    AUC         = {}
    dist = np.arange(0,200.5)
    
    colors    = {'kaf':'b', 'naf':'g'}
    for i,indx in enumerate(models[sort_i_count]):
        listOfList.append([])
        for chan in excDict:
            variables = model_sets[indx]['variables'][chan]
            if chan not in AUC: AUC[chan]=[]
            if chan == 'kaf': a4=1
            elif chan == 'naf': a4=1-variables[1]
            a5=variables[1]
            a6=variables[2]
            a7=variables[3]
            auc_func = np.power(10,variables[0])*(a4 + a5/(1 + np.exp((dist-a6)/a7) ))
            AUC[chan].append( np.trapz(auc_func) )
            for j,val in enumerate(model_sets[indx]['variables'][chan]):
                if i == 0: excDict[chan][j] = []
                excDict[chan][j].append( val )
                listOfList[i].append(val)
    
    
    for chan in excDict:
        fig,ax = plt.subplots(1,1, figsize=(3,3))
        ax.plot(np.array(cs_proportion)[cs_pos_models], np.array(AUC[chan])[cs_pos_models], 
                            'o', ms=10, mew=1, mec='w', c=colors[chan])
        plot_linreg(ax, np.array(cs_proportion)[cs_pos_models], np.array(AUC[chan])[cs_pos_models])
        fig.savefig('Figures/lincorr_csProp-AUC{}-cond{}_ramping{}.png'.format(condition,chan,tag), dpi=300, transparent=True)
        
    plt.show()
    #plot_AUC(model_sets)
    
    
    '''
    for chan in excDict:
        fchan, achan = plt.subplots( max(2,len(excDict[chan])), 1, figsize=(18, 3*max(2,len(excDict[chan]))) )
        achan[0].set_title(chan)
        for l,val in excDict[chan].items():
            achan[l].plot(range(54), val, 'k')
    ''' 
    #cluster_kmean( np.array(listOfList) )
    #pca( np.array(listOfList) )



def plot_linreg(ax,X,Y):
    from scipy.stats import linregress
    x       = np.array(X)
    out     = linregress(X,Y)
    y       = out[0]*x + out[1]
    ax.plot(x,y)
    ax.set_title('r={:.1f}, p={:.1e}'.format(out[2],out[3]), pad=0, fontsize=20 )
    xmin, xmax, ymin, ymax = ax.axis()
    add_box2ax(ax, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    ax.axis('off')
    #ax.invert_xaxis()
    #ax.invert_yaxis()
    
    print(out)



def df_proportion_cs(df, ax=None):
    
    if not ax: fig, ax = plt.subplots(1,1)
        
    df_negative = df.loc[df['cs'] == 0]
    df_pos      = df.loc[df['cs'] == 1]
    
    n1,b1,p1 = ax.hist(    df['condition'], 
                range=(-0.5,3.5), 
                bins=4, 
                color='lightgrey', 
                histtype='step',
                lw=4, 
                alpha=1.0
                )
    n2,b2,p2 = ax.hist(    df_pos['condition'], 
                range=(-0.5,3.5), 
                bins=4, 
                color='r', 
                histtype='stepfilled', 
                alpha=1.0
                )
    
    percent = n2/(n1)*100
    
    for i,x in enumerate([2,3]):
        ax.text(x, n2[x]+100, '{0:.0f}%'.format(percent[x]))
    
    ax.set_xticks([0,1,2,3])
    ax.set_xticklabels(['ctrl','ACh','DA','DA+ACh'])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)



def df_model_cs_ratio(df, ax=None):
        
    df_negative = df.loc[df['cs'] == 0]
    df_pos      = df.loc[df['cs'] == 1]
    
    nbins = df.model.max()
    count = 3
    
    df_pos_cond = df_pos.loc[df_pos['condition'] == count]
    
    # sort 
    models  = np.array(sorted(df.model.unique()))
    cs_proportion  = []
    
    for m in models:
        cs_proportion.append( len(df_pos_cond.loc[df_pos_cond['model']==m]) )
    sort_i_count    = np.argsort(cs_proportion)
    
    nbins           = len(models)
    theList         = []
    for i,indx in enumerate(sort_i_count):
        theList += [i]*cs_proportion[indx]
        
    if not ax: fig, ax = plt.subplots(1,1, figsize=(18,9))
    
    n2,b2,p2 = ax.hist( theList, 
                range=(-0.5,nbins-0.5),  
                bins=nbins, 
                color='r', 
                histtype='bar', 
                alpha=1.0,
                ec='black'
                )
    
    # sortera modeller på mest cs 
    ax.set_xticks(range(len(n2)))
    ax.set_xticklabels( models[sort_i_count], fontsize='small' )
    
    # correlate with model parameters ========================================================================================
    path            = '../Cholinergic_modulation'
    with open('{}/D1_71bestFit_updRheob.pkl'.format(path), 'rb') as f:
        model_sets  = pickle.load(f, encoding="latin1")
    
    # rheobase ----------------------------------------------------
    excList = []
    fmain, amain = plt.subplots( 1, 1, figsize=(18,9) )
    for i,indx in enumerate(models[sort_i_count]):
        #print(indx)
        excList.append(model_sets[indx]['rheobase'])    #model_sets[cell_index]['variables']
    #plt.figure(figsize=(18,9))
    ax.fill_between([0,54], [100,100], color='lightgrey')
    ax.plot([0,54], [50,50], 'w')
    amain.plot(range(54), np.array(excList)/max(excList), 'r', zorder=1e6)
    ax.plot(range(54), np.array(excList)/max(excList)*(100), 'k', lw=5)
    amain.set_xticks([])
    
    # model variables ----------------------------------------------------
    excDict = {'naf':{}, 'kaf':{}}#, 'kas':{}, 'kir':{}, 'sk':{}} #, 'can':{}, 'c32':{}, 'c33':{}}
    listOfList = []
    AUC = {}
    dist = np.arange(0,200.5)
    fkaf,akaf = plt.subplots(1,1)
    fnaf,anaf = plt.subplots(1,1)
    achandist = {'kaf':akaf, 'naf':anaf}
    lw=3
    for A in achandist.values():
        A.plot([0,200],[0,0],'k',lw=lw, dash_capstyle='butt')
        A.plot([100,100],[-0.02,0.02],'w',lw=lw, zorder=5e7)
        A.set_xlim([0,200])
    
    
    colors    = {'kaf':'b', 'naf':'g'}
    for i,indx in enumerate(models[sort_i_count]):
        listOfList.append([])
        for chan in excDict:
            variables = model_sets[indx]['variables'][chan]
            if chan not in AUC: AUC[chan]=[]
            if chan == 'kaf': a4=1
            elif chan == 'naf': a4=1-variables[1]
            a5=variables[1]
            a6=variables[2]
            a7=variables[3]
            auc_func = np.power(10,variables[0])*(a4 + a5/(1 + np.exp((dist-a6)/a7) ))
            if i == 0: achandist[chan].fill_between(dist, auc_func, color=colors[chan])
            achandist[chan].plot( dist, auc_func, 'k')
            AUC[chan].append( np.trapz(auc_func) )
            for j,val in enumerate(model_sets[indx]['variables'][chan]):
                if i == 0: excDict[chan][j] = []
                excDict[chan][j].append( val )
                listOfList[i].append(val)
     
    ax.plot(range(54), np.array(AUC['kaf'])/max(AUC['kaf'])*(100), color=colors['kaf'], lw=5, zorder=2e6)
    ax.plot(range(54), np.array(AUC['naf'])/max(AUC['naf'])*(100), color=colors['naf'], lw=5, zorder=2e6)
    
    # plot
    
    for chan in excDict:
        fchan, achan = plt.subplots( max(2,len(excDict[chan])), 1, figsize=(18, 3*max(2,len(excDict[chan]))) )
        achan[0].set_title(chan)
        for l,val in excDict[chan].items():
            achan[l].plot(range(54), val, 'k')
        
        fmain.savefig('Figures/cs_proportion-all-model-var-{}.png'.format(chan), dpi=300, transparent=True)
    
    ax.axis('off'); akaf.axis('off'); anaf.axis('off')
    fig.savefig('Figures/cs_proportion-model-cond{}.png'.format(['ctrl','ACh','DA','ACh+DA'][count]), dpi=300, transparent=True)
    fmain.savefig('Figures/cs_proportion-model-cond{}_rheobase.png'.format(['ctrl','ACh','DA','ACh+DA'][count]), dpi=300, transparent=True)
    fnaf.savefig('Figures/chan_dist_models_naf.png', dpi=300, transparent=True)
    fkaf.savefig('Figures/chan_dist_models_kaf.png', dpi=300, transparent=True)
    
    #cluster_kmean( np.array(listOfList) )
    #pca( np.array(listOfList) )

def cluster_kmean(X):    
    ''' atempt to find multidimentional relationships.
    Not well tested and nothing interesting came out.'''
    
    # cluster (Kmean) -----------------------
    from sklearn.cluster import KMeans
    # Number of clusters
    nclust = 10
    kmeans = KMeans(n_clusters=nclust)
    # Fitting the input data
    kmeans = kmeans.fit(X)
    # Getting the cluster labels
    labels = kmeans.labels_
    # Centroid values
    centroids = kmeans.cluster_centers_

def pca(X):     
    ''' atempt to find multidimentional relationships.
    Not well tested and nothing interesting came out.'''
    fk, ak = plt.subplots(1,1,figsize=(18,9))
    ak.plot(range(len(n2)), labels, 'ow', mew=2, mec='k')
    ak.set_title('kmin')
    ak.set_xticks(range(len(n2)))
    ak.set_xticklabels( models[sort_i_count], fontsize='small' )
    
    # 
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    x = StandardScaler().fit_transform(X)
    pca = PCA()
    x_pca = pca.fit_transform(x)
    x_pca = pd.DataFrame(x_pca)
    x_pca['target'] = labels
    labels = []
    for i in range(X.shape[1]):
        labels.append('PC{}'.format(i+1))
    x_pca.columns = labels+['target']
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('Principal Component 1') 
    ax.set_ylabel('Principal Component 2') 
    targets = range(nclust)
    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'] #['r', 'g', 'b']
    for target, color in zip(targets,colors):
     indicesToKeep = x_pca['target'] == target
     ax.scatter(x_pca.loc[indicesToKeep, 'PC1']
     , x_pca.loc[indicesToKeep, 'PC2']
     , c = color
     , s = 50)
    ax.legend(targets)
    ax.grid()
    
            
def extract_spikes(load=False):
    
    if load:
        # load and serialize
        with open('Results/extracted_spikes.json', 'r') as handle:
            spikes = json.load(handle)
        return np.array(spikes)
        
    # (number of trials in data; up to 15 first spikes)
    nTrials = 54*5*(1+2+10+20)
    spike_array = np.zeros((nTrials,15))
    trialID = 0
    for c in conditions:
        file_string = 'Results/inVivo_ramping_{}_model*.json'.format(c) 
        files = glob.glob(file_string)
        for f in files:
            with open(f, 'r') as handle:
                data = json.load(handle)
            # levels in dict: fid_ach -> fid_da -> condition -> bg
            # e.g. trace = data[fid_ach][fid_da][condition][bg]
            for key,achlev in data.items():
                if key == 'time': continue
                for dalev in achlev.values():
                    for d in dalev[c].values():
                        
                        # get spikes
                        spike_train = use.getSpikedata_x_y(data['time'],d)
                        spike_array[trialID, :len(spike_train)] = spike_train
                        trialID += 1
    
    # remove columns at end with 0 spikes
    ind = -1; flag = False
    while np.sum(spike_array[:,ind]) == 0:
        ind -= 1; flag=1
        print(ind)
    
    if flag: 
        spike_array = spike_array[:,:ind]
    
    # serialize and save
    with open('Results/extracted_spikes.json', 'wt') as handle:
        json.dump(spike_array.tolist(), handle, indent=4)
        
    return spike_array       

         
def extract_data_into_frame(load=False, tag='', path='Results/'):
    
    if load:
        # load and serialize
        return pd.read_pickle('{}cs_extracted_into_df{}.pkl'.format(path,tag))
        
    # (number of trials in data; up to 15 first spikes)
    columns = ['model', 'condition', 'bgid', '#spikes', 'cs', 'fACh', 'fDA']
    nTrials = 54*5*(1+2+10+20)
    df      = pd.DataFrame(columns=columns, index=range(nTrials)).T
    trialID = 0
    for ci,c in enumerate(conditions):
        file_string = 'Results/inVivo{}_ramping_{}_model*.json'.format(tag,c) 
        files = glob.glob(file_string)
        print(c, len(files))
        for f in files:
            with open(f, 'r') as handle:
                data = json.load(handle)
            model = int(f.split('model')[1].split('.')[0])
            # levels in dict: fid_ach -> fid_da -> condition -> bg
            # e.g. trace = data[fid_ach][fid_da][condition][bg]
            for key,achlev in data.items():
                if key == 'time': continue
                for i,dalev in enumerate(achlev.values()):
                    for bg,d in enumerate(dalev[c].values()):
                        # get spikes
                        spike_train = use.getSpikedata_x_y(data['time'],d)
                        # cs?
                        b = f4a.check_sliding_average_lowSampleData(d, threshold=-37)
                        if b: w = 1
                        else: w = 0
                        df[trialID] = [model, c, bg, len(spike_train), w, int(key), i]
                        trialID += 1
    
    # serialize and save
    df = df.T
    df.to_pickle('{}cs_extracted_into_df{}.pkl'.format(path,tag))
        
    return df     

def plot_spike_hist(spike_array, ax=None):
    '''
    plots histogram of spikes under all conditions
    '''
    # TODO
    # update conditions to be more reliable, 
    #   e.g. work with different number of trails than the hard coded ones.
    if not ax: fig,ax = plt.subplots(1,1)
    
    index = 0; prevind = 0
    for i,ii in enumerate([1,2,10,20]):
        index = index + nbg*ii*nmodels
        
        irev = int(20/ii)
        # get only lines corresponding to one condition
        sub_array   = np.repeat(spike_array[prevind:index,:]-1000, irev)
        prevind     = index
        
        # plot hist without zeros [sub_array != 0]
        ax.hist(    sub_array[sub_array>0], 
                    bins=15, 
                    range= (minT,maxT),
                    density=0, 
                    color=colors[conditions[i]],
                    alpha=0.8,
                    zorder=1000+ontop[i]
                    )
        ax.hist(    sub_array[sub_array>0], 
                    bins=15, 
                    range= (minT,maxT),
                    density=0, 
                    histtype='step',
                    color=colors[conditions[i]],
                    zorder=ontop[i]
                    )
    ax.axis('off')
    
def plot_spike_rast(spike_array, ax=None):
    '''
    raster plot...
    '''
    if not ax: fig,ax = plt.subplots(1,1)
    
    gap = 500
    
    index = 0; prevind = 0
    for i,ii in enumerate([1,2,10,20]):
        index = index + nbg*ii*nmodels
        
        # get only lines corresponding to one condition
        sub_array   = spike_array[prevind:index,:]-1000
        
        # fix line values (y)
        y = np.ones(sub_array.shape)
        for l in range(y.shape[0]):
            y[l,:] = prevind + gap*i + l
            
        plt.plot(sub_array, y, '.', ms=2, c=colors[conditions[i]])
        
        prevind = index
        
        ax.plot([minT,maxT], [index+gap*(i+0.5), index+gap*(i+0.5)], '--k', lw=0.5)  
    ax.set_ylim([-0.5*gap,index+gap*(i+0.5)])

# -------------------------------------------------------------------------------------------
# 3rd

def add_grey_bg(ax,tmin=200,ymax=40,ymin=-70,yl1=-30,yl2=-50,fillcol='#E0E0E0'):
    
    ax.fill_between([tmin,maxT], [ymax,ymax], y2=[ymin,ymin], color=fillcol)
    ax.plot([tmin,maxT], [yl1,yl1], '-w')
    ax.plot([tmin,maxT], [yl2,yl2], '-w')
    ax.plot([300,300], [ymin,ymax], '-w')
    ax.plot([400,400], [ymin,ymax], '-w')


def add_box2ax(ax, xmin=0, ymin=0, xmax=100, ymax=100, lw=2, ls='-', c='k'):
    
    ax.plot([xmin,xmax],[ymin,ymin], c=c, ls=ls, lw=lw)
    ax.plot([xmin,xmax],[ymax,ymax], c=c, ls=ls, lw=lw)
    ax.plot([xmin,xmin],[ymin,ymax], c=c, ls=ls, lw=lw)
    ax.plot([xmax,xmax],[ymin,ymax], c=c, ls=ls, lw=lw)


def get_ax():
    fig     = plt.figure(figsize=(6,4))
    r       = 6
    grid    = plt.GridSpec(r, r, hspace=0.1, wspace=0)
    a1      = fig.add_subplot(grid[:,:-1])
    a2      = fig.add_subplot(grid[:,-1])
    return fig, [a1,a2]

def get_candidates(save2File=False):
    
    cand = [{'mod':2, 'fach':0, 'fda':4, 'bg':4},
            {'mod':2, 'fach':0, 'fda':8, 'bg':4},
            {'mod':2, 'fach':0, 'fda':9, 'bg':4},
            {'mod':2, 'fach':2, 'fda':1, 'bg':4},
            {'mod':2, 'fach':2, 'fda':3, 'bg':4},
            {'mod':2, 'fach':2, 'fda':8, 'bg':2},
            {'mod':5, 'fach':0, 'fda':0, 'bg':0},
            {'mod':5, 'fach':0, 'fda':0, 'bg':4},
            {'mod':5, 'fach':0, 'fda':1, 'bg':4},
            {'mod':5, 'fach':0, 'fda':3, 'bg':4},
            {'mod':5, 'fach':0, 'fda':4, 'bg':4},
            {'mod':5, 'fach':0, 'fda':5, 'bg':4},
            {'mod':5, 'fach':0, 'fda':6, 'bg':2},
            {'mod':5, 'fach':0, 'fda':7, 'bg':4},
            {'mod':5, 'fach':0, 'fda':8, 'bg':4},
            {'mod':5, 'fach':0, 'fda':9, 'bg':4},
            {'mod':5, 'fach':2, 'fda':0, 'bg':4},
            {'mod':5, 'fach':2, 'fda':1, 'bg':4},
            {'mod':5, 'fach':2, 'fda':2, 'bg':1},
            {'mod':5, 'fach':2, 'fda':3, 'bg':4},
            {'mod':5, 'fach':2, 'fda':4, 'bg':2},
            {'mod':5, 'fach':2, 'fda':7, 'bg':2},
            {'mod':5, 'fach':2, 'fda':7, 'bg':4},
            {'mod':5, 'fach':2, 'fda':8, 'bg':2},
            {'mod':5, 'fach':2, 'fda':8, 'bg':4},
            {'mod':10, 'fach':0, 'fda':0, 'bg':0},
            {'mod':10, 'fach':0, 'fda':0, 'bg':1},
            {'mod':10, 'fach':0, 'fda':0, 'bg':2},
            {'mod':10, 'fach':0, 'fda':0, 'bg':3},
            {'mod':10, 'fach':0, 'fda':0, 'bg':4},
            {'mod':10, 'fach':0, 'fda':1, 'bg':0},
            {'mod':10, 'fach':0, 'fda':1, 'bg':3},
            {'mod':10, 'fach':0, 'fda':1, 'bg':4},
            {'mod':10, 'fach':0, 'fda':2, 'bg':0},
            {'mod':10, 'fach':0, 'fda':2, 'bg':3},
            {'mod':10, 'fach':0, 'fda':2, 'bg':4},
            {'mod':10, 'fach':0, 'fda':3, 'bg':0},
            {'mod':10, 'fach':0, 'fda':3, 'bg':1},
            {'mod':10, 'fach':0, 'fda':3, 'bg':2},
            {'mod':10, 'fach':0, 'fda':3, 'bg':3},
            {'mod':10, 'fach':0, 'fda':3, 'bg':4},
            {'mod':10, 'fach':0, 'fda':4, 'bg':0},
            {'mod':10, 'fach':0, 'fda':4, 'bg':1},
            {'mod':10, 'fach':0, 'fda':4, 'bg':2},
            {'mod':10, 'fach':0, 'fda':4, 'bg':3},
            {'mod':10, 'fach':0, 'fda':4, 'bg':4},
            {'mod':10, 'fach':0, 'fda':5, 'bg':0},
            {'mod':10, 'fach':0, 'fda':5, 'bg':3},
            {'mod':10, 'fach':0, 'fda':5, 'bg':4},
            {'mod':10, 'fach':0, 'fda':6, 'bg':1},
            {'mod':10, 'fach':0, 'fda':6, 'bg':2},
            {'mod':10, 'fach':0, 'fda':6, 'bg':3},
            {'mod':10, 'fach':0, 'fda':6, 'bg':4},
            {'mod':10, 'fach':0, 'fda':7, 'bg':0},
            {'mod':10, 'fach':0, 'fda':7, 'bg':1},
            {'mod':10, 'fach':0, 'fda':7, 'bg':2},
            {'mod':10, 'fach':0, 'fda':7, 'bg':3},
            {'mod':10, 'fach':0, 'fda':7, 'bg':4},
            {'mod':10, 'fach':0, 'fda':8, 'bg':0},
            {'mod':10, 'fach':0, 'fda':8, 'bg':1},
            {'mod':10, 'fach':0, 'fda':8, 'bg':2},
            {'mod':10, 'fach':0, 'fda':8, 'bg':3},
            {'mod':10, 'fach':0, 'fda':8, 'bg':4},
            {'mod':10, 'fach':0, 'fda':9, 'bg':0},
            {'mod':10, 'fach':0, 'fda':9, 'bg':1},
            {'mod':10, 'fach':0, 'fda':9, 'bg':3},
            {'mod':10, 'fach':0, 'fda':9, 'bg':4},
            {'mod':10, 'fach':2, 'fda':0, 'bg':0},
            {'mod':10, 'fach':2, 'fda':0, 'bg':1},
            {'mod':10, 'fach':2, 'fda':0, 'bg':2},
            {'mod':10, 'fach':2, 'fda':0, 'bg':3},
            {'mod':10, 'fach':2, 'fda':0, 'bg':4},
            {'mod':10, 'fach':2, 'fda':1, 'bg':0},
            {'mod':10, 'fach':2, 'fda':1, 'bg':1},
            {'mod':10, 'fach':2, 'fda':1, 'bg':2},
            {'mod':10, 'fach':2, 'fda':1, 'bg':3},
            {'mod':10, 'fach':2, 'fda':1, 'bg':4},
            {'mod':10, 'fach':2, 'fda':2, 'bg':1},
            {'mod':10, 'fach':2, 'fda':2, 'bg':2},
            {'mod':10, 'fach':2, 'fda':2, 'bg':3},
            {'mod':10, 'fach':2, 'fda':2, 'bg':4},
            {'mod':10, 'fach':2, 'fda':3, 'bg':0},
            {'mod':10, 'fach':2, 'fda':3, 'bg':1},
            {'mod':10, 'fach':2, 'fda':3, 'bg':2},
            {'mod':10, 'fach':2, 'fda':3, 'bg':3},
            {'mod':10, 'fach':2, 'fda':3, 'bg':4},
            {'mod':10, 'fach':2, 'fda':4, 'bg':1},
            {'mod':10, 'fach':2, 'fda':4, 'bg':2},
            {'mod':10, 'fach':2, 'fda':4, 'bg':3},
            {'mod':10, 'fach':2, 'fda':4, 'bg':4},
            {'mod':10, 'fach':2, 'fda':5, 'bg':0},
            {'mod':10, 'fach':2, 'fda':5, 'bg':1},
            {'mod':10, 'fach':2, 'fda':5, 'bg':2},
            {'mod':10, 'fach':2, 'fda':5, 'bg':3},
            {'mod':10, 'fach':2, 'fda':5, 'bg':4},
            {'mod':10, 'fach':2, 'fda':6, 'bg':1},
            {'mod':10, 'fach':2, 'fda':6, 'bg':2},
            {'mod':10, 'fach':2, 'fda':6, 'bg':3},
            {'mod':10, 'fach':2, 'fda':6, 'bg':4},
            {'mod':10, 'fach':2, 'fda':7, 'bg':0},
            {'mod':10, 'fach':2, 'fda':7, 'bg':1},
            {'mod':10, 'fach':2, 'fda':7, 'bg':2},
            {'mod':10, 'fach':2, 'fda':7, 'bg':3},
            {'mod':10, 'fach':2, 'fda':7, 'bg':4},
            {'mod':10, 'fach':2, 'fda':8, 'bg':0},
            {'mod':10, 'fach':2, 'fda':8, 'bg':1},
            {'mod':10, 'fach':2, 'fda':8, 'bg':2},
            {'mod':10, 'fach':2, 'fda':8, 'bg':3},
            {'mod':10, 'fach':2, 'fda':8, 'bg':4},
            {'mod':10, 'fach':2, 'fda':9, 'bg':0},
            {'mod':10, 'fach':2, 'fda':9, 'bg':2},
            {'mod':10, 'fach':2, 'fda':9, 'bg':3},
            {'mod':10, 'fach':2, 'fda':9, 'bg':4},
            {'mod':15, 'fach':0, 'fda':0, 'bg':4},
            {'mod':15, 'fach':0, 'fda':1, 'bg':4},
            {'mod':15, 'fach':0, 'fda':3, 'bg':4},
            {'mod':15, 'fach':0, 'fda':4, 'bg':3},
            {'mod':15, 'fach':0, 'fda':4, 'bg':4},
            {'mod':15, 'fach':0, 'fda':5, 'bg':4},
            {'mod':15, 'fach':0, 'fda':6, 'bg':4},
            {'mod':15, 'fach':0, 'fda':7, 'bg':4},
            {'mod':15, 'fach':0, 'fda':8, 'bg':4},
            {'mod':15, 'fach':2, 'fda':0, 'bg':4},
            {'mod':15, 'fach':2, 'fda':1, 'bg':0},
            {'mod':15, 'fach':2, 'fda':1, 'bg':4},
            {'mod':15, 'fach':2, 'fda':2, 'bg':4},
            {'mod':15, 'fach':2, 'fda':3, 'bg':2},
            {'mod':15, 'fach':2, 'fda':3, 'bg':4},
            {'mod':15, 'fach':2, 'fda':4, 'bg':3},
            {'mod':15, 'fach':2, 'fda':4, 'bg':4},
            {'mod':15, 'fach':2, 'fda':5, 'bg':0},
            {'mod':15, 'fach':2, 'fda':5, 'bg':4},
            {'mod':15, 'fach':2, 'fda':6, 'bg':4},
            {'mod':15, 'fach':2, 'fda':7, 'bg':4},
            {'mod':15, 'fach':2, 'fda':8, 'bg':2},
            {'mod':15, 'fach':2, 'fda':8, 'bg':3},
            {'mod':15, 'fach':2, 'fda':8, 'bg':4},
            {'mod':15, 'fach':2, 'fda':9, 'bg':4},
            {'mod':27, 'fach':0, 'fda':2, 'bg':4},
            {'mod':27, 'fach':0, 'fda':3, 'bg':4},
            {'mod':27, 'fach':0, 'fda':4, 'bg':4},
            {'mod':27, 'fach':0, 'fda':6, 'bg':4},
            {'mod':27, 'fach':2, 'fda':3, 'bg':4},
            {'mod':27, 'fach':2, 'fda':4, 'bg':2},
            {'mod':27, 'fach':2, 'fda':4, 'bg':4},
            {'mod':27, 'fach':2, 'fda':5, 'bg':4},
            {'mod':27, 'fach':2, 'fda':7, 'bg':4},
            {'mod':28, 'fach':2, 'fda':3, 'bg':2},
            {'mod':28, 'fach':2, 'fda':7, 'bg':2},
            {'mod':33, 'fach':0, 'fda':3, 'bg':4},
            {'mod':33, 'fach':0, 'fda':7, 'bg':4},
            {'mod':39, 'fach':0, 'fda':0, 'bg':4},
            {'mod':39, 'fach':2, 'fda':1, 'bg':4},
            {'mod':39, 'fach':2, 'fda':5, 'bg':4},
            {'mod':39, 'fach':2, 'fda':6, 'bg':4},
            {'mod':47, 'fach':0, 'fda':6, 'bg':4},
            {'mod':47, 'fach':0, 'fda':8, 'bg':4},
            {'mod':47, 'fach':2, 'fda':3, 'bg':4},
            {'mod':47, 'fach':2, 'fda':7, 'bg':4},
            {'mod':47, 'fach':2, 'fda':8, 'bg':4}
            ] 
    
    if save2File:
        with open('Results/example_candidates_cs_ramping.json', 'wt') as handle:
            json.dump(cand, handle, indent=4)
    
    return cand


def get_candidates_full(save2file=False):
    
    cand = [{'mod':2, 'fach':0, 'fda':4, 'bg':4},
            {'mod':2, 'fach':0, 'fda':8, 'bg':4},
            {'mod':2, 'fach':0, 'fda':9, 'bg':4},
            {'mod':2, 'fach':2, 'fda':1, 'bg':4},
            {'mod':2, 'fach':2, 'fda':3, 'bg':4},
            {'mod':2, 'fach':2, 'fda':8, 'bg':2},
            {'mod':5, 'fach':0, 'fda':0, 'bg':0},
            {'mod':5, 'fach':0, 'fda':0, 'bg':4},
            {'mod':5, 'fach':0, 'fda':1, 'bg':4},
            {'mod':5, 'fach':0, 'fda':3, 'bg':4},
            {'mod':5, 'fach':0, 'fda':4, 'bg':4},
            {'mod':5, 'fach':0, 'fda':5, 'bg':4},
            {'mod':5, 'fach':0, 'fda':6, 'bg':2},
            {'mod':5, 'fach':0, 'fda':7, 'bg':4},
            {'mod':5, 'fach':0, 'fda':8, 'bg':4},
            {'mod':5, 'fach':0, 'fda':9, 'bg':4},
            {'mod':5, 'fach':2, 'fda':0, 'bg':4},
            {'mod':5, 'fach':2, 'fda':1, 'bg':4},
            {'mod':5, 'fach':2, 'fda':2, 'bg':1},
            {'mod':5, 'fach':2, 'fda':3, 'bg':4},
            {'mod':5, 'fach':2, 'fda':4, 'bg':2},
            {'mod':5, 'fach':2, 'fda':7, 'bg':2},
            {'mod':5, 'fach':2, 'fda':7, 'bg':4},
            {'mod':5, 'fach':2, 'fda':8, 'bg':2},
            {'mod':5, 'fach':2, 'fda':8, 'bg':4},
            {'mod':10, 'fach':0, 'fda':0, 'bg':0},
            {'mod':10, 'fach':0, 'fda':0, 'bg':1},
            {'mod':10, 'fach':0, 'fda':0, 'bg':2},
            {'mod':10, 'fach':0, 'fda':0, 'bg':3},
            {'mod':10, 'fach':0, 'fda':0, 'bg':4},
            {'mod':10, 'fach':0, 'fda':1, 'bg':0},
            {'mod':10, 'fach':0, 'fda':1, 'bg':3},
            {'mod':10, 'fach':0, 'fda':1, 'bg':4},
            {'mod':10, 'fach':0, 'fda':2, 'bg':0},
            {'mod':10, 'fach':0, 'fda':2, 'bg':3},
            {'mod':10, 'fach':0, 'fda':2, 'bg':4},
            {'mod':10, 'fach':0, 'fda':3, 'bg':0},
            {'mod':10, 'fach':0, 'fda':3, 'bg':1},
            {'mod':10, 'fach':0, 'fda':3, 'bg':2},
            {'mod':10, 'fach':0, 'fda':3, 'bg':3},
            {'mod':10, 'fach':0, 'fda':3, 'bg':4},
            {'mod':10, 'fach':0, 'fda':4, 'bg':0},
            {'mod':10, 'fach':0, 'fda':4, 'bg':1},
            {'mod':10, 'fach':0, 'fda':4, 'bg':2},
            {'mod':10, 'fach':0, 'fda':4, 'bg':3},
            {'mod':10, 'fach':0, 'fda':4, 'bg':4},
            {'mod':10, 'fach':0, 'fda':5, 'bg':0},
            {'mod':10, 'fach':0, 'fda':5, 'bg':3},
            {'mod':10, 'fach':0, 'fda':5, 'bg':4},
            {'mod':10, 'fach':0, 'fda':6, 'bg':1},
            {'mod':10, 'fach':0, 'fda':6, 'bg':2},
            {'mod':10, 'fach':0, 'fda':6, 'bg':3},
            {'mod':10, 'fach':0, 'fda':6, 'bg':4},
            {'mod':10, 'fach':0, 'fda':7, 'bg':0},
            {'mod':10, 'fach':0, 'fda':7, 'bg':1},
            {'mod':10, 'fach':0, 'fda':7, 'bg':2},
            {'mod':10, 'fach':0, 'fda':7, 'bg':3},
            {'mod':10, 'fach':0, 'fda':7, 'bg':4},
            {'mod':10, 'fach':0, 'fda':8, 'bg':0},
            {'mod':10, 'fach':0, 'fda':8, 'bg':1},
            {'mod':10, 'fach':0, 'fda':8, 'bg':2},
            {'mod':10, 'fach':0, 'fda':8, 'bg':3},
            {'mod':10, 'fach':0, 'fda':8, 'bg':4},
            {'mod':10, 'fach':0, 'fda':9, 'bg':0},
            {'mod':10, 'fach':0, 'fda':9, 'bg':1},
            {'mod':10, 'fach':0, 'fda':9, 'bg':3},
            {'mod':10, 'fach':0, 'fda':9, 'bg':4},
            {'mod':10, 'fach':2, 'fda':0, 'bg':0},
            {'mod':10, 'fach':2, 'fda':0, 'bg':1},
            {'mod':10, 'fach':2, 'fda':0, 'bg':2},
            {'mod':10, 'fach':2, 'fda':0, 'bg':3},
            {'mod':10, 'fach':2, 'fda':0, 'bg':4},
            {'mod':10, 'fach':2, 'fda':1, 'bg':0},
            {'mod':10, 'fach':2, 'fda':1, 'bg':1},
            {'mod':10, 'fach':2, 'fda':1, 'bg':2},
            {'mod':10, 'fach':2, 'fda':1, 'bg':3},
            {'mod':10, 'fach':2, 'fda':1, 'bg':4},
            {'mod':10, 'fach':2, 'fda':2, 'bg':1},
            {'mod':10, 'fach':2, 'fda':2, 'bg':2},
            {'mod':10, 'fach':2, 'fda':2, 'bg':3},
            {'mod':10, 'fach':2, 'fda':2, 'bg':4},
            {'mod':10, 'fach':2, 'fda':3, 'bg':0},
            {'mod':10, 'fach':2, 'fda':3, 'bg':1},
            {'mod':10, 'fach':2, 'fda':3, 'bg':2},
            {'mod':10, 'fach':2, 'fda':3, 'bg':3},
            {'mod':10, 'fach':2, 'fda':3, 'bg':4},
            {'mod':10, 'fach':2, 'fda':4, 'bg':1},
            {'mod':10, 'fach':2, 'fda':4, 'bg':2},
            {'mod':10, 'fach':2, 'fda':4, 'bg':3},
            {'mod':10, 'fach':2, 'fda':4, 'bg':4},
            {'mod':10, 'fach':2, 'fda':5, 'bg':0},
            {'mod':10, 'fach':2, 'fda':5, 'bg':1},
            {'mod':10, 'fach':2, 'fda':5, 'bg':2},
            {'mod':10, 'fach':2, 'fda':5, 'bg':3},
            {'mod':10, 'fach':2, 'fda':5, 'bg':4},
            {'mod':10, 'fach':2, 'fda':6, 'bg':1},
            {'mod':10, 'fach':2, 'fda':6, 'bg':2},
            {'mod':10, 'fach':2, 'fda':6, 'bg':3},
            {'mod':10, 'fach':2, 'fda':6, 'bg':4},
            {'mod':10, 'fach':2, 'fda':7, 'bg':0},
            {'mod':10, 'fach':2, 'fda':7, 'bg':1},
            {'mod':10, 'fach':2, 'fda':7, 'bg':2},
            {'mod':10, 'fach':2, 'fda':7, 'bg':3},
            {'mod':10, 'fach':2, 'fda':7, 'bg':4},
            {'mod':10, 'fach':2, 'fda':8, 'bg':0},
            {'mod':10, 'fach':2, 'fda':8, 'bg':1},
            {'mod':10, 'fach':2, 'fda':8, 'bg':2},
            {'mod':10, 'fach':2, 'fda':8, 'bg':3},
            {'mod':10, 'fach':2, 'fda':8, 'bg':4},
            {'mod':10, 'fach':2, 'fda':9, 'bg':0},
            {'mod':10, 'fach':2, 'fda':9, 'bg':2},
            {'mod':10, 'fach':2, 'fda':9, 'bg':3},
            {'mod':10, 'fach':2, 'fda':9, 'bg':4},
            {'mod':11, 'fach':2, 'fda':0, 'bg':4},
            {'mod':15, 'fach':0, 'fda':0, 'bg':4},
            {'mod':15, 'fach':0, 'fda':1, 'bg':4},
            {'mod':15, 'fach':0, 'fda':3, 'bg':4},
            {'mod':15, 'fach':0, 'fda':4, 'bg':3},
            {'mod':15, 'fach':0, 'fda':4, 'bg':4},
            {'mod':15, 'fach':0, 'fda':5, 'bg':4},
            {'mod':15, 'fach':0, 'fda':6, 'bg':4},
            {'mod':15, 'fach':0, 'fda':7, 'bg':4},
            {'mod':15, 'fach':0, 'fda':8, 'bg':4},
            {'mod':15, 'fach':2, 'fda':0, 'bg':4},
            {'mod':15, 'fach':2, 'fda':1, 'bg':0},
            {'mod':15, 'fach':2, 'fda':1, 'bg':4},
            {'mod':15, 'fach':2, 'fda':2, 'bg':4},
            {'mod':15, 'fach':2, 'fda':3, 'bg':2},
            {'mod':15, 'fach':2, 'fda':3, 'bg':4},
            {'mod':15, 'fach':2, 'fda':4, 'bg':3},
            {'mod':15, 'fach':2, 'fda':4, 'bg':4},
            {'mod':15, 'fach':2, 'fda':5, 'bg':0},
            {'mod':15, 'fach':2, 'fda':5, 'bg':4},
            {'mod':15, 'fach':2, 'fda':6, 'bg':4},
            {'mod':15, 'fach':2, 'fda':7, 'bg':4},
            {'mod':15, 'fach':2, 'fda':8, 'bg':2},
            {'mod':15, 'fach':2, 'fda':8, 'bg':3},
            {'mod':15, 'fach':2, 'fda':8, 'bg':4},
            {'mod':15, 'fach':2, 'fda':9, 'bg':4},
            {'mod':17, 'fach':0, 'fda':3, 'bg':4},
            {'mod':17, 'fach':0, 'fda':4, 'bg':4},
            {'mod':17, 'fach':0, 'fda':6, 'bg':4},
            {'mod':17, 'fach':0, 'fda':7, 'bg':4},
            {'mod':17, 'fach':0, 'fda':8, 'bg':4},
            {'mod':17, 'fach':2, 'fda':0, 'bg':4},
            {'mod':17, 'fach':2, 'fda':1, 'bg':4},
            {'mod':17, 'fach':2, 'fda':2, 'bg':4},
            {'mod':17, 'fach':2, 'fda':5, 'bg':4},
            {'mod':17, 'fach':2, 'fda':8, 'bg':4},
            {'mod':27, 'fach':0, 'fda':2, 'bg':4},
            {'mod':27, 'fach':0, 'fda':3, 'bg':4},
            {'mod':27, 'fach':0, 'fda':4, 'bg':4},
            {'mod':27, 'fach':0, 'fda':6, 'bg':4},
            {'mod':27, 'fach':2, 'fda':3, 'bg':4},
            {'mod':27, 'fach':2, 'fda':4, 'bg':2},
            {'mod':27, 'fach':2, 'fda':4, 'bg':4},
            {'mod':27, 'fach':2, 'fda':5, 'bg':4},
            {'mod':27, 'fach':2, 'fda':7, 'bg':4},
            {'mod':28, 'fach':2, 'fda':3, 'bg':2},
            {'mod':28, 'fach':2, 'fda':7, 'bg':2},
            {'mod':33, 'fach':0, 'fda':3, 'bg':4},
            {'mod':33, 'fach':0, 'fda':7, 'bg':4},
            {'mod':39, 'fach':0, 'fda':0, 'bg':4},
            {'mod':39, 'fach':2, 'fda':1, 'bg':4},
            {'mod':39, 'fach':2, 'fda':5, 'bg':4},
            {'mod':39, 'fach':2, 'fda':6, 'bg':4},
            {'mod':47, 'fach':0, 'fda':6, 'bg':4},
            {'mod':47, 'fach':0, 'fda':8, 'bg':4},
            {'mod':47, 'fach':2, 'fda':3, 'bg':4},
            {'mod':47, 'fach':2, 'fda':7, 'bg':4},
            {'mod':47, 'fach':2, 'fda':8, 'bg':4},
            {'mod':56, 'fach':0, 'fda':0, 'bg':4},
            {'mod':56, 'fach':0, 'fda':1, 'bg':4},
            {'mod':56, 'fach':0, 'fda':2, 'bg':4},
            {'mod':56, 'fach':0, 'fda':3, 'bg':1},
            {'mod':56, 'fach':0, 'fda':3, 'bg':4},
            {'mod':56, 'fach':0, 'fda':4, 'bg':1},
            {'mod':56, 'fach':0, 'fda':4, 'bg':4},
            {'mod':56, 'fach':0, 'fda':5, 'bg':4},
            {'mod':56, 'fach':0, 'fda':6, 'bg':4},
            {'mod':56, 'fach':0, 'fda':7, 'bg':1},
            {'mod':56, 'fach':0, 'fda':7, 'bg':4},
            {'mod':56, 'fach':0, 'fda':8, 'bg':4},
            {'mod':56, 'fach':2, 'fda':0, 'bg':1},
            {'mod':56, 'fach':2, 'fda':0, 'bg':2},
            {'mod':56, 'fach':2, 'fda':0, 'bg':3},
            {'mod':56, 'fach':2, 'fda':0, 'bg':4},
            {'mod':56, 'fach':2, 'fda':1, 'bg':2},
            {'mod':56, 'fach':2, 'fda':1, 'bg':4},
            {'mod':56, 'fach':2, 'fda':2, 'bg':2},
            {'mod':56, 'fach':2, 'fda':2, 'bg':4},
            {'mod':56, 'fach':2, 'fda':3, 'bg':1},
            {'mod':56, 'fach':2, 'fda':3, 'bg':2},
            {'mod':56, 'fach':2, 'fda':3, 'bg':4},
            {'mod':56, 'fach':2, 'fda':4, 'bg':1},
            {'mod':56, 'fach':2, 'fda':4, 'bg':3},
            {'mod':56, 'fach':2, 'fda':4, 'bg':4},
            {'mod':56, 'fach':2, 'fda':5, 'bg':2},
            {'mod':56, 'fach':2, 'fda':5, 'bg':4},
            {'mod':56, 'fach':2, 'fda':6, 'bg':1},
            {'mod':56, 'fach':2, 'fda':6, 'bg':4},
            {'mod':56, 'fach':2, 'fda':7, 'bg':2},
            {'mod':56, 'fach':2, 'fda':7, 'bg':3},
            {'mod':56, 'fach':2, 'fda':7, 'bg':4},
            {'mod':56, 'fach':2, 'fda':8, 'bg':0},
            {'mod':56, 'fach':2, 'fda':8, 'bg':1},
            {'mod':56, 'fach':2, 'fda':8, 'bg':2},
            {'mod':56, 'fach':2, 'fda':8, 'bg':3},
            {'mod':56, 'fach':2, 'fda':8, 'bg':4},
            {'mod':56, 'fach':2, 'fda':9, 'bg':0},
            {'mod':56, 'fach':2, 'fda':9, 'bg':1},
            {'mod':56, 'fach':2, 'fda':9, 'bg':2},
            {'mod':56, 'fach':2, 'fda':9, 'bg':4},
            {'mod':63, 'fach':0, 'fda':5, 'bg':4},
            {'mod':69, 'fach':0, 'fda':8, 'bg':4},
            {'mod':69, 'fach':2, 'fda':3, 'bg':4},
            {'mod':69, 'fach':2, 'fda':5, 'bg':4},
            {'mod':69, 'fach':2, 'fda':7, 'bg':2}
            ]
    
    if save2file:
        with open('../all_cs_candidates_ramping.json', 'wt') as handle:
            json.dump(cand, handle, indent=4)
    
    return cand
        

run()
