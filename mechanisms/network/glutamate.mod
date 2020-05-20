COMMENT
Updated Exp2Syn synapse with Mg-blocked nmda channel.

Defaul values of parameters (time constants etc) set to match synaptic channels in 
striatal medium spiny neurons (Du et al., 2017; Chapman et al., 2003; Ding et al., 2008).

Robert . Lindroos @ ki . se

original comment:
________________
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT






NEURON {
	POINT_PROCESS glutamate
	RANGE tau1_ampa, tau2_ampa, tau1_nmda, tau2_nmda
	RANGE erev, g, i
	RANGE i_ampa, i_nmda, g_ampa, g_nmda, ratio, I, G, mg, q, block, alpha, beta
	RANGE ampa_scale_factor, nmda_scale_factor
    RANGE damod, maxModNMDA,max2NMDA,maxModAMPA,max2AMPA,l1NMDA,l2NMDA,l1AMPA,l2AMPA
	
	NONSPECIFIC_CURRENT i
	USEION cal WRITE ical VALENCE 2
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	erev        = 0.0       (mV)
	
	tau1_ampa   = 1.9       (ms)
    tau2_ampa   = 4.8       (ms)  : tau2 > tau1
    tau1_nmda   = 5.52      (ms)  : Chapman et al 2003; table 1, adult rat (rise time, rt = 12.13. rt ~= 2.197*tau (wiki;rise time) -> tau = 12.13 / 2.197 ~= 5.52
    tau2_nmda   = 231       (ms)  : Chapman et al 2003 (table 1; adult)
    
    ratio       = 1         (1)   : both components give same maximal amplitude of current
    mg          = 1         (mM)
    alpha       = 0.062
    beta        = 3.57
    q           = 2               : approx room temp -> 
    
    nmda_scale_factor = 1
    ampa_scale_factor = 1
    
    ca_ratio_ampa = 0.005
    ca_ratio_nmda = 0.1
    
    maxModNMDA  = 1
    max2NMDA    = 1
    maxModAMPA  = 1
    max2AMPA    = 1
    damod       = 0
    l1NMDA      = 0
    l2NMDA      = 0
    l1AMPA      = 0
    l2AMPA      = 0
}


ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor_nmda
	factor_ampa
	i_ampa
	i_nmda
	g_ampa
	g_nmda
	block
	I
	G
	ical (nA)
}


STATE {
	A (uS)
	B (uS)
	C (uS)
	D (uS)
}



INITIAL {
	LOCAL tp
	if (tau1_nmda/tau2_nmda > .9999) {
		tau1_nmda = .9999*tau2_nmda
	}
	if (tau1_ampa/tau2_ampa > .9999) {
		tau1_ampa = .9999*tau2_ampa
	}
	
	: NMDA
	A           = 0
	B           = 0
	tp          = (tau1_nmda*tau2_nmda)/(tau2_nmda - tau1_nmda) * log(tau2_nmda/tau1_nmda)
	factor_nmda = -exp(-tp/tau1_nmda) + exp(-tp/tau2_nmda)
	factor_nmda = 1/factor_nmda
	
	: AMPA
	C           = 0
	D           = 0
	tp          = (tau1_ampa*tau2_ampa)/(tau2_ampa - tau1_ampa) * log(tau2_ampa/tau1_ampa)
	factor_ampa = -exp(-tp/tau1_ampa) + exp(-tp/tau2_ampa)
	factor_ampa = 1/factor_ampa
}




BREAKPOINT {
	SOLVE state METHOD cnexp
	
	: NMDA
	g_nmda = (B - A) * modulation(maxModNMDA,max2NMDA,l1NMDA,l2NMDA)
	block  = MgBlock()
	i_nmda = g_nmda * (v - erev) * block * nmda_scale_factor
	
	: AMPA
	g_ampa = (D - C) * modulation(maxModAMPA,max2AMPA,l1AMPA,l2AMPA)
	i_ampa = g_ampa * (v - erev) * ampa_scale_factor
	
	: total current
	G = g_ampa + g_nmda
	I = i_ampa + i_nmda
	
	: splitting in ca and non ca currents
	ical = i_ampa*ca_ratio_ampa  + i_nmda*ca_ratio_nmda
    i = i_ampa*(1-ca_ratio_ampa) + i_nmda*(1-ca_ratio_nmda)
}



DERIVATIVE state {
	A' = -A/tau1_nmda*q
	B' = -B/tau2_nmda*q
	C' = -C/tau1_ampa*q
	D' = -D/tau2_ampa*q
}



NET_RECEIVE(weight (uS)) {
	A = A + weight*factor_nmda
	B = B + weight*factor_nmda
	C = C + weight*factor_ampa*ratio
	D = D + weight*factor_ampa*ratio
}



FUNCTION MgBlock() {
    
    MgBlock = 1 / (1 + mg * exp(-alpha * v) / beta )
    
}

FUNCTION modulation(m1,m2,l1,l2) {
    : returns modulation factor
    
    modulation = 1 + damod * ( (m1-1)*l1 + (m2-1)*l2 )
    if (modulation < 0) {
        modulation = 0
    } 
}




