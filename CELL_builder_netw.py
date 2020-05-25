#
'''
The MSN class defining the cell
    updated based on MSN_builder in base of repo to handle network models (specifically
    different paramfile structure)

the ffactor is used to rescale leak and capacitance when adding spines. 
It is the reversed factors used in Du et al. 2017.
'''

from neuron import h
import numpy as np
import json
from math import exp

import logging
logger = logging.getLogger(__name__)

FLOAT_FORMAT = '%.17g'

# ======================= the MSN class ==================================================
class CELL:
    def __init__(self,  params=None,
                        morphology=None,
                        mechanisms=None,
                        variables=None,
                        replace_axon=True,
                        N=40.0,
                        ffactor=None
                        ):
        Import = h.Import3d_SWC_read()
        Import.input(morphology)
        imprt = h.Import3d_GUI(Import, 0)
        imprt.instantiate(None)
        h.define_shape()
        
        self._set_nsegs(N=N)
        self._create_sectionlists(replace_axon=replace_axon)
        self._read_param_file(params)
        
        # initialize soma as start point of distance function
        #h.distance(sec=self.soma)
        # Define origin of distance function
        h.distance(0, 0.5, sec=self.soma)
        
        # global params---hard coded for now. TODO find better solution
        h.celsius = self.par['global']['celsius']['value']
        self.v_init = self.par['global']['v_init']['value']
        
        if mechanisms:
            with open(mechanisms) as file:
                mechLists = json.load(file)
             
        
        # set channels and parameters
        for region, seclist in zip(['somatic', 'axonal', 'basal'],[self.somalist, self.axonlist, self.dendlist]):
            
            # if no mechanism file is passed assumes cadyn objects for spn. 
            if mechanisms:              cl = mechLists[region]
            elif region == 'axonal':    cl = self.channel_lists[region]
            else:                       cl = self.channel_lists[region]+["cadyn_ms", "caldyn_ms"]
            
            for sec in seclist:
                for mech in cl:
                    sec.insert(mech)
                
                if not mechanisms or ('all' in mechLists and 'pas' in mechLists['all']):
                    sec.insert('pas')
                    sec.e_pas = self.par['section'][region]['e_pas']['value']
                    if ffactor:
                        sec.g_pas = self.par['section'][region]['g_pas']['value'] / ffactor
                    else:
                        sec.g_pas = self.par['section'][region]['g_pas']['value']
                    
                # passive params (section params)---hard coded for now. TODO find better solution
                sec.Ra = self.par['section'][region]['Ra']['value']
                if ffactor:
                    sec.cm = self.par['section'][region]['cm']['value'] / ffactor
                else:
                    sec.cm = self.par['section'][region]['cm']['value']
                sec.ena = self.par['section'][region]['ena']['value']
                sec.ek = self.par['section'][region]['ek']['value']
                
                # channels (range params)
                for rp in self.par['range'][region].keys():
                    
                    value = self.par['range'][region][rp]['value']
                    
                    if self.par['range'][region][rp]['dist_type'] in ['exp','distance']:
                        function = self.par['range'][region][rp]['dist']
                        # get dist
                        for seg in sec:
                            distance = h.distance(1, seg.x, sec=sec)
                            val = eval( function.format(distance=distance, value=value) )
                            
                            cmd = 'seg.{} = {}'.format(rp, val)
                            exec(cmd)
                    else:
                        val = value 
                        for seg in sec:
                            cmd = 'seg.{} = {}'.format(rp, val)
                            exec(cmd)
                            
        
    def _read_param_file(self, params):
        par             = { "global":   {},
                            "section":  {"axonal":{}, "somatic":{}, "basal":{}},
                            "range":    {"axonal":{}, "somatic":{}, "basal":{}} }
        channel_lists   = { "axonal":[], "somatic":[], "basal":[] }
        
        with open(params) as file:
            # OBS. hard coded to use first set of modulation parameters
            # TODO update?
            paramList = json.load(file)[0]
        
        for p in paramList:
            # extract values etc
            name    = p["param_name"]
            value   = p["value"]
            
            ptype   = p["type"]
            if ptype == 'global':
                par["global"][name] = {'value':value}
            else:
                sl      = p["sectionlist"]
                dt      = p["dist_type"]
                par[ptype][sl][name] = {'value':value, 'dist_type':dt}
                if ptype == 'range':
                    mech = p["mech"]
                    channel_lists[sl].append(mech)
                    if 'dist' in p:
                        par[ptype][sl][name]['dist'] = p['dist']
                        
                
        self.par = par
        self.channel_lists = channel_lists
            
        
    def _create_sectionlists(self, replace_axon=False):
        # soma
        self.nsomasec = 0
        self.somalist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('soma') >= 0:
                self.somalist.append(sec=sec)
                if self.nsomasec == 0:
                    self.soma = sec
                self.nsomasec += 1
        # dendrite
        self.dendlist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('dend') >= 0:
                self.dendlist.append(sec=sec)
        # axon
        self.axonlist = h.SectionList()
        if replace_axon:
            self._create_AIS()
        else:
            axon=[]
            for sec in h.allsec():
                if sec.name().find('axon') >= 0:
                    self.axonlist.append(sec=sec)
        # all
        self.allsecnames = []
        self.allseclist  = h.SectionList()
        for sec in h.allsec():
            self.allsecnames.append(sec.name())
            self.allseclist.append(sec=sec)
        
        
    
    
    def _set_nsegs(self, N=40.0):
        """ def seg/sec """
        
        for sec in h.allsec():
            if 'dend' in sec.name():
                #print('inne', sec.name())
                sec.nseg = 2*int(sec.L/N)+1
            else:
                sec.nseg = 1
            
    
    def _create_AIS(self):
        """Replica of "Replace axon" in: 
            https://bluepyopt.readthedocs.io/en/latest/_modules/bluepyopt/ephys/morphologies.html#Morphology
        """
        
        temp = []
        for sec in h.allsec():
            if sec.name().find('axon') >= 0:
                temp.append(sec)
        
        # specify diameter based on blu
        if len(temp) == 0:
            ais_diams = [1, 1]
        elif len(temp) == 1:
            ais_diams = [temp[0].diam, temp[0].diam]
        else:
            ais_diams = [temp[0].diam, temp[0].diam]
            # Define origin of distance function
            h.distance(0, 0.5, sec=self.soma)
            
            for section in h.allsec():
                if section.name().find('axon') >= 0:
                    # If distance to soma is larger than 60, store diameter
                    if h.distance(1, 0.5, sec=section) > 60:
                        ais_diams[1] = section.diam
                        break
        
        # delete old axon
        for section in temp:
            h.delete_section(sec=section)
        
        # Create new axon sections
        a0 = h.Section(name='axon[0]')
        a1 = h.Section(name='axon[1]')
        
        # populate axonlist
        for sec in h.allsec():
            if sec.name().find('axon') >= 0:
                self.axonlist.append(sec=sec)
        
        # connect axon sections to soma and eachother
        a0.connect(self.soma)
        a1.connect(a0)
        
        # set axon params
        for index, section in enumerate([a0,a1]):
            section.nseg = 1
            section.L = 30
            section.diam = ais_diams[index]
        
        # this line is needed to prevent garbage collection of axon 
        self.axon = [a0,a1]
        
        logger.debug('Replace axon with AIS') 
        
        

    
                    
