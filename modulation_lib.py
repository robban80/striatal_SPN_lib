


    

class DA():
    ''' 
    class for setting DA modulation of cell. 
    arguments:
        cell = cell object to modulate
        mod_dict = dict with chan:factor pairs
    additional arguments carries default behavior of modulation.
    
    The 'naf', 'kas', 'kir', 'cal12','cal13','can' and 'car channels are reported to be modulated by DA in dspn. Aditionally the car channel is modulated in ispn.
    
    By default the AIS is not modulated, since this was found to increase chance of positive modulation in dspn.
        ---OBS. if uniform modulation is wanted this has to be stated!!!!
    
    Modulation has in some channel(s?) been reported to be a shift in gating. This is not
        implemented for DA (but see ACh class for example).
    
    Modulation can be set individually for gaba, glut and intrinsic.
        by default all are used.
    
    The _reset_mod method can be used to turn of modulation.
    '''
    
    def __init__(self,  cell, mod_dict,                                    
                        modulation='noAxon', 
                        intrinsic_mod=1,
                        gaba_mod=1,
                        glut_mod=1,
                        play=[],
                        syn_dict={},
                        dt=0.025
                        ):
        
        self.cell = cell
        self.mod_dict = mod_dict
        self.syn_dict = syn_dict
        self.play = play
        self.dt = dt
        
        if modulation == 'uniform':
            self.compartments = [cell.dendlist, cell.somalist, cell.axonlist]
        elif modulation == 'noAxon':
            self.compartments = [cell.dendlist, cell.somalist]
        else: raise Exception('Error, "{}" modulation, not permitted\n\tuse "uniform" or "noAxon"'.format(modulation))
        
        
        self._set_modulation(   intrinsic_mod,
                                gaba_mod,
                                glut_mod
                                )
        
    def _set_modulation(self,
                        intrinsic_mod,
                        gaba_mod,
                        glut_mod):
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    if intrinsic_mod:   self._update_conductance(seg)
                    if gaba_mod:        self._set_gaba(seg)
                    if glut_mod:        self._set_glut(seg)
    
    
    def _update_conductance(self, seg, reset=0):
        for mech in seg:
            if mech.name() in self.mod_dict:
                # set shift
                mech.damod = 1
                mech.maxMod = self.mod_dict[mech.name()]
                if reset: 
                    mech.level = 0
                elif len(self.play) and mech.name() in self.play:
                    self.play[mech.name()].play(mech._ref_level, self.dt)
                else:
                    mech.level = 1
                            
    def _set_gaba(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'gaba' in syn.hname():
                syn.damod = 1
                syn.maxMod = self.syn_dict['GABA']
                if reset:
                    syn.level = 0 
                elif len(self.play) and 'gaba' in self.play:
                    self.play['gaba'].play(syn._ref_level, self.dt)
                else:
                    syn.level = 1
    
    def _set_glut(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'glut' in syn.hname():
                syn.damod = 1
                syn.maxModNMDA = self.syn_dict['NMDA']
                syn.maxModAMPA = self.syn_dict['AMPA']
                if reset:
                    syn.l1AMPA = 0; syn.l1NMDA = 0
                elif len(self.play) and 'glut' in self.play:
                    self.play['glut'].play(syn._ref_l1NMDA, self.dt)
                    self.play['glut'].play(syn._ref_l1AMPA, self.dt)
                else:
                    syn.l1AMPA = 1; syn.l1NMDA = 1
    
    def _reset_mod(self):
        if len(self.play):
            for trans in self.play.values():
                trans.play_remove()
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    self._update_conductance(seg, reset=1)
                    self._set_glut(seg, reset=1)
                    self._set_gaba(seg, reset=1)
      

class ACh():
    ''' 
    class for setting ACh modulation of cell. 
    arguments:
        cell = cell object to modulate
        mod_dict = dict with chan:factor pairs
    additional arguments carries default behavior of modulation.
    
    The 'naf', 'kir','cal12','cal13','can' and 'im' channels are reported to be modulated by ACh.
    Additionally "kaf" is (by default) shifted 20 mv more negative
    
    Modulation can be set individually for gaba, glut, intrinsic and kaf shift.
        by default all are used.
    '''
    
    def __init__(self,  cell, mod_dict,                                    
                        modulation='uniform', 
                        intrinsic_mod=1,
                        shift_kaf=1,
                        gaba_mod=1,
                        glut_mod=1,
                        mv_shift_kaf=20,
                        syn_dict={},
                        play=[],
                        dt=0.025
                        ):
        
        self.cell = cell
        self.mod_dict = mod_dict
        self.mv_shift_kaf = mv_shift_kaf
        self.syn_dict = syn_dict
        self.play = play
        self.dt = dt
        
        if modulation == 'uniform':
            self.compartments = [cell.dendlist, cell.somalist, cell.axonlist]
        elif modulation == 'noAxon':
            self.compartments = [cell.dendlist, cell.somalist]
        else: raise Exception('Error, "{}" modulation, not permitted\n\tuse "uniform" or "noAxon"'.format(modulation))
        
        self._set_modulation(   intrinsic_mod,
                                shift_kaf,
                                gaba_mod,
                                glut_mod
                                )
        
    def _set_modulation(self,
                        intrinsic_mod,
                        shift_kaf,
                        gaba_mod,
                        glut_mod):
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    if intrinsic_mod:   self._update_conductance(seg)
                    if shift_kaf:       self._shift_kaf(seg)
                    if gaba_mod:        self._set_gaba(seg)
                    if glut_mod:        self._set_glut(seg)
        
    def _update_conductance(self, seg, reset=0):
        for mech in seg:
            if mech.name() in self.mod_dict:
                # set shift
                mech.damod = 1
                mech.max2 = self.mod_dict[mech.name()]
                if reset:
                    mech.lev2 = 0
                if len(self.play) and mech.name() in self.play:
                    self.play[mech.name()].play(mech._ref_lev2, self.dt)
                else:
                    mech.lev2 = 1
    
    def _shift_kaf(self, seg, reset=0):
        for mech in seg:
            if mech.name() == 'kaf':
                # set shift
                if reset:
                    mech.modShift = 0
                elif len(self.play) and 'kaf' in self.play:
                    self.play['kaf'].play(mech._ref_modShift, self.dt)
                else:
                    mech.modShift = self.mv_shift_kaf
                            
    def _set_gaba(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'gaba' in syn.hname():
                syn.damod = 1
                syn.max2 = self.syn_dict['GABA']
                if reset:
                    syn.lev2 = 0
                elif len(self.play) and 'gaba' in self.play:
                    self.play['gaba'].play(syn._ref_lev2, self.dt)
                else:
                    syn.lev2 = 1
    
    def _set_glut(self, seg, reset=0):
        for syn in seg.point_processes():
            if 'glut' in syn.hname():
                syn.damod = 1
                syn.max2NMDA = self.syn_dict['NMDA']
                syn.max2AMPA = self.syn_dict['AMPA']
                if reset:
                    syn.l2AMPA = 0
                    syn.l2NMDA = 0
                elif len(self.play) and 'glut' in self.play:
                    self.play['glut'].play(syn._ref_l2NMDA, self.dt)
                    self.play['glut'].play(syn._ref_l2AMPA, self.dt)
                else:
                    syn.l2AMPA = 1
                    syn.l2NMDA = 1
                            
    
    def _reset_mod(self):
        if len(self.play):
            for trans in self.play.values():
                trans.play_remove()
        for comp in self.compartments:
            for sec in comp:
                for seg in sec:
                    self._shift_kaf(seg, reset=1)
                    self._update_conductance(seg, reset=1)
                    self._set_glut(seg, reset=1)
                    self._set_gaba(seg, reset=1)
    
