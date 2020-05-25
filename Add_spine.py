
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).

#from neuron import h

class Spine():
    """
    Spine class. Creates a spine with neck and head.
    based on Mattioni and Le Novere, (2013).
    https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=150284&file=/TimeScales-master/neuronControl/spine.py#tabs-2
    
    The spine will be connected to section at location x
    
    To move a spine use method move_spine(new_parent, x)
    """
    
    def __init__(self, h,               \
                       parent_section,  \
                       x,               \
                       h0=0,            \
                       n0=0,            \
                       neck_L=1.0,      \
                       neck_dia=0.1,    \
                       head_L=0.5,      \
                       head_dia=0.5,    \
                       Ra=150.0,        \
                       passive=False,   \
                       Cm=1.0           ):
        """ Create a spine with geometry given by the arguments
            neck cytoplasmic resistivity set to 795 ohm-cm giving a resistance of 500 Mohm
            (using default geometry) following Dorman et al., 2018 -> Harnett et al., 2012
            *Synaptic amplification by dendritic spines enhances input cooperativity*
            """
        
        self.h          =   h
        self.parent     =   parent_section  # the parent section connected to the neck
        self.x          =   x
        self.passive    =   passive
        self.neck       =   self.create_neck(neck_L=neck_L, neck_dia=neck_dia, Cm=Cm, Ra=785)
        self.head       =   self.create_head(head_L=head_L, head_dia=head_dia, Cm=Cm, Ra=Ra)
        self.stim       =   None    # attribute for saving spike apparatus
        self.input2head =   h0      # flag for checking if synaptic input is located in head
        self.input2neck =   n0      # flag for checking if synaptic input is located in neck
        
        # connect head to neck
        self.head.connect(self.neck(1),0)
        
        # connect neck to parent section
        self.neck.connect(parent_section(x),0)
        
        # calculate somatic distance of parent seg
        # TODO: how do we know what the reference point is? 
        #       Should be cell->soma, and seems to be in test (since same dist obtained as in Alex models)
        #       -but can we know for sure that this is always the case?
        #       https://www.neuron.yale.edu/neuron/static/py_doc/modelspec/programmatic/topology/geometry.html#distance
        self.setDistance()
        
    
    def create_neck(self, neck_L=1.0, neck_dia=0.1, Ra=785.0, Cm=1.0):
        """ Create a spine neck"""
        
        # TODO should the sec_name be unique?
        sec_name        =   'spine_neck'
        neck            =   self.h.Section(name=sec_name)
        neck.nseg       =   1
        neck.L          =   neck_L 
        neck.diam       =   neck_dia
        neck.Ra         =   Ra 
        neck.cm         =   Cm
        
        if not self.passive:
            for mech in [   'cat32',    \
                            'cat33'     ]:
                neck.insert('{}_ms'.format(mech))
        
            neck(0.5).pbar_cat32_ms = 1e-7
            neck(0.5).pbar_cat33_ms = 1e-8
        
        neck.insert('pas')
        neck.g_pas      =   1.25e-5
        neck.e_pas      =   -70 
        neck.insert('caldyn_ms')
        
        return neck
        
        
        
    def create_head(self, head_L=0.5, head_dia=0.5, Ra=150.0, Cm=1.0):
        """Create the head of the spine and populate it with channels"""
        
        # TODO should the sec_name be unique?
        sec_name        =   'spine_head'
        head            =   self.h.Section(name=sec_name)
        
        head.nseg       =   1
        head.L          =   head_L
        head.diam       =   head_dia
        head.Ra         =   Ra
        head.cm         =   1.0
        
        if not self.passive:
            for mech in [   'kir',      \
                            'cat32',    \
                            'cat33',    \
                            'car',      \
                            'cal12',    \
                            'cal13',    \
                            'cadyn'     ]:
                head.insert('{}_ms'.format(mech))
            
            head(0.5).pbar_cat32_ms = 1e-7
            head(0.5).pbar_cat33_ms = 1e-8
            head(0.5).pbar_car_ms   = 1e-8
            head(0.5).pbar_cal12_ms = 1e-7
            head(0.5).pbar_cal13_ms = 1e-8
            head(0.5).gbar_kir_ms   = 1e-7
            head.ek                 = -90
        
        head.insert('pas')
        head.g_pas      =   1.25e-5
        head.e_pas      =   -70 
        head.insert('caldyn_ms')
        
        
        return head
    
    def setDistance(self):
        ''' set distance of spine to soma
            '''
        dist = self.h.distance(self.x, sec=self.parent)
        #if dist > 90 and dist < 120 and 'dSPN' in self.parent.name():
        #    print(self.parent.name(), self.x, dist)
        self.somaticDist = dist
        
    def move_spine(self, new_parent, x):
        ''' move spine from one section to another
        '''
        self.h.disconnect(sec=self.neck)
        self.neck.connect(new_parent(x),0)
        self.parent = new_parent
        self.x      = x
        
