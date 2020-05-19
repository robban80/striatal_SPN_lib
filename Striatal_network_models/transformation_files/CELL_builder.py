#
'''
The MSN class defining the cell
    updated based on MSN_builder in base of repo to handle network models (specifically
    different paramfile structure)
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
                        section=None
                        ):
        Import = h.Import3d_SWC_read()
        Import.input(morphology)
        imprt = h.Import3d_GUI(Import, 0)
        imprt.instantiate(None)
        h.define_shape()
        
        self._set_nsegs()
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
                    sec.g_pas = self.par['section'][region]['g_pas']['value']
                    
                # passive params (section params)---hard coded for now. TODO find better solution
                sec.Ra = self.par['section'][region]['Ra']['value']
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
        
        
    
    
    def _set_nsegs(self):
        """ def seg/sec """
        
        for sec in h.allsec():
            sec.nseg = 2*int(sec.L/40.0)+1
            
    
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
        
        

                           
   

class Spine():
    """
    Spine class. Create a spine with neck and head.
    inspired by Mattioni and Le Novere, (2013).
    https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=150284&file=/TimeScales-master/neuronControl/spine.py#tabs-2
    
    if a parent section is passed as argument, the spine will be connected to that section (using the default orientation)
        to connect a spine at a later stage use the connect_spine(parent) method, where parent is a section.
        
    Since default orientation is used, to place multiple spines spread over one section, first use function split_section(sec) in common_functions.
        split_section(sec) splits one section into multiple sections with retained total surface area (and close conductances).
    """
    #TODO: 
    # -test spine (one and multiple): rheobase etc retained?
    # -run in vivo with spines
    
    def __init__(self, id,              \
                       neck_L=1.0,      \
                       neck_dia=0.1,    \
                       head_L=0.5,      \
                       head_dia=0.5,    \
                       Ra=150.0,        \
                       Cm=1.0,          \
                       parent=None      ):
        """ Create a spine with geometry given by the arguments"""
        
        self.id         =   id
        self.neck       =   self.create_neck(neck_L=neck_L, neck_dia=neck_dia, Ra=Ra, Cm=Cm)
        self.head       =   self.create_head(head_L=head_L, head_dia=head_dia, Ra=Ra, Cm=Cm)
        self.parent     =   parent  # the parent section connected to the neck
        self.stim       =   None    # attribute for saving spike apparatus (netStim, synapse and
        
        self.connect_head2neck(parent)
        
    
    def create_neck(self, neck_L=1.0, neck_dia=0.1, Ra=150.0, Cm=1.0):
        """ Create a spine neck"""
        
        sec_name        =   'spine_%d_neck' % (self.id)
        neck            =   h.Section(name=sec_name)
        neck.nseg       =   1
        neck.L          =   neck_L 
        neck.diam       =   neck_dia
        neck.Ra         =   Ra 
        neck.cm         =   Cm
        
        for mech in [   'pas',      \
                        'cav32',    \
                        'cav33',    \
                        'caldyn'     ]:
            neck.insert(mech)
        
        neck(0.5).pbar_cav33 = 1e-7
        neck(0.5).pbar_cav33 = 1e-8
        
        neck.g_pas      =   1.25e-5
        neck.e_pas      =   -70 
        
        return neck
        
        
        
    def create_head(self, head_L=0.5, head_dia=0.5, Ra=150.0, Cm=1.0):
        """Create the head of the spine and populate it with channels"""
        
        sec_name        =   'spine_%d_head' % (self.id)
        head            =   h.Section(name=sec_name)
        
        head.nseg       =   1
        head.L          =   head_L
        head.diam       =   head_dia
        head.Ra         =   Ra
        head.cm         =   1.0
        
        for mech in [   'pas',      \
                        'kir',      \
                        'cav32',    \
                        'cav33',    \
                        'car',      \
                        'cal12',    \
                        'cal13',    \
                        'cadyn',    \
                        'caldyn'    ]:
            head.insert(mech)
        
        head(0.5).pbar_cav32 = 1e-7
        head(0.5).pbar_cav33 = 1e-8
        head(0.5).pbar_car   = 1e-8
        head(0.5).pbar_cal12 = 1e-7
        head(0.5).pbar_cal13 = 1e-8
        head(0.5).gbar_kir   = 1e-7
        
        head.g_pas      =   1.25e-5
        head.e_pas      =   -70 
        
        return head
        
    
    def connect_head2neck(self, parent=None):
        ''' connect spine head to neck and if parent is not None connect spine to parent 
        connection is hardcoded to use default connection orientation (0 end connected to 1 end on parent)
        To use other orientation, first create and then connect spine, 
            using the connect spine method with updated arguments
        '''
        self.head.connect(self.neck(1),0)
        if not parent is None:
            self.neck.connect(parent(1),0) 
    
    def connect_spine(self, parent, x=1, end=0):
        ''' connect spine to parent sec.
        '''
        self.neck.connect(parent(x),end)
        self.parent = parent
    
    def move_spine(self, new_parent):
        ''' move spine from one section to another (using default orientation)
        '''
        h.disconnect(sec=self.neck)
        self.neck.connect(new_parent(1),0)
        self.parent = new_parent
        
           
    
    
                    
