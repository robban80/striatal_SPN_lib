TITLE Inwardly rectifying potassium current

COMMENT
neuromodulation is added as functions:
    
    modulation = 1 + damod*(maxMod-1)*level

where:
    
    damod  [0]: is a switch for turning modulation on or off {1/0}
    maxMod [1]: is the maximum modulation for this specific channel (read from the param file)
                e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}
    level  [0]: is an additional parameter for scaling modulation. 
                Can be used simulate non static modulation by gradually changing the value from 0 to 1 {0-1}

[] == default values
{} == ranges
    
ENDCOMMENT

NEURON {
    THREADSAFE
    SUFFIX kir
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik
    RANGE damod, maxMod, level, max2, lev2
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 3
    damod = 0
    maxMod = 1
    level = 0
    max2 = 1
    lev2 = 0
} 

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    gk (S/cm2)
    minf
    mtau (ms)
}

STATE { m }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gbar*m*modulation()
    ik = gk*(v-ek)
}

DERIVATIVE states {
    rates()
    m' = (minf-m)/mtau*q
}

INITIAL {
    rates()
    m = minf
}

PROCEDURE rates() {
    LOCAL alpha, beta, sum
    UNITSOFF
    minf = 1/(1+exp((v-(-102))/13))
    alpha = 0.1*exp((v-(-60))/(-14))
    beta = 0.27/(1+exp((v-(-31))/(-23)))
    sum = alpha+beta
    mtau = 1/sum
    UNITSON
}

FUNCTION modulation() {
    : returns modulation factor
    
    modulation = 1 + damod * ( (maxMod-1)*level + (max2-1)*lev2 ) 
    if (modulation < 0) {
        modulation = 0
    }  
}

COMMENT

Original data by Steephen (2009), rat, room temp.

Genesis implementation by Kai Du <kai.du@ki.se>, MScell v9.5.

NEURON implementation by Alexander Kozlov <akozlov@csc.kth.se>, smooth
fit of mtau.

ENDCOMMENT
