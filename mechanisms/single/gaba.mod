COMMENT
Updated Exp2Syn synapse including a modulation function.

Defaul values of parameters (time constants etc) set to match synaptic channels in 
striatal medium spiny neurons, following Wolf et al., 2007 -> Galaretta 1997

    modulation = 1 + damod*( (maxMod-1)*level + (max2-1)*lev2 )

two substrates can modulate independently, e.g. DA and ACh

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
	POINT_PROCESS gaba
	RANGE tau1, tau2
	RANGE erev, g, i, q
    RANGE damod, maxMod, level, max2, lev2
	
	NONSPECIFIC_CURRENT i
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	erev    = -60.0     (mV)
	tau1    =   0.5     (ms)    : wolf et al., 2007 -> Galaretta 1997
    tau2    =   7.5     (ms)    : wolf et al., 2007 -> Galaretta 1997
    q       =   2       : approx room temp 
    
    damod       = 0
    maxMod      = 1
    max2        = 1
    level       = 0
    lev2        = 0
}


ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
}


STATE {
	A (uS)
	B (uS)
}


INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A       = 0
	B       = 0
	tp      = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor  = -exp(-tp/tau1) + exp(-tp/tau2)
	factor  = 1/factor
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	
	g = (B - A) * modulation(maxMod,max2,level,lev2)
	i = g * (v - erev)
}


DERIVATIVE state {
	A' = -A/tau1*q
	B' = -B/tau2*q
}


NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}


FUNCTION modulation(m1,m2,l1,l2) {
    : calculates modulation factor
    
    modulation = 1 + damod * ( (m1-1)*l1 + (m2-1)*l2 )
    if (modulation < 0) {
        modulation = 0
    } 
}




