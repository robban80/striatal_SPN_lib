TITLE Fast A-type potassium current (Kv4.2)

COMMENT
neuromodulation is added as in two ways: 
1. A shift parameter is added to the infinity curves that can be used to shift the voltage
    dependece of the channel. By deafult this does not also shift the time constant.
2. A function that scales the conductance of the channel (voltage independent).


functions:
    
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
    SUFFIX kaf
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik
    RANGE damod, maxMod, level, max2, lev2, modShift
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 2
    damod = 0
    maxMod = 1
    level = 0
    max2 = 1
    lev2 = 0
    modShift = 0
} 

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    gk (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gbar*m*m*h *modulation()
    ik = gk*(v-ek)
}

DERIVATIVE states {
    rates()
    m' = (minf-m)/mtau*q
    h' = (hinf-h)/htau*q
}

INITIAL {
    rates()
    m = minf
    h = hinf
}

PROCEDURE rates() {
    LOCAL alpha, beta, sum
    UNITSOFF
    alpha = 1.5/(1+exp((v-4+modShift)/(-17)))
    beta = 0.6/(1+exp((v-10+modShift)/9))
    sum = alpha+beta
    minf = alpha/sum
    mtau = 1/sum
    : mtau = 1/( 1.5/(1+exp((v-4)/(-17))) + 0.6/(1+exp((v-10)/9)) ) : don't shift tau
    
    alpha = 0.105/(1+exp((v-(-121)+modShift)/22))
    beta = 0.065/(1+exp((v-(-55)+modShift)/(-11)))
    sum = alpha+beta
    hinf = alpha/sum
    htau = 1/sum
    : htau = 1/( 0.105/(1+exp((v-(-121))/22)) + 0.065/(1+exp((v-(-55))/(-11))) ) : don't shift tau
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

Original data by Tkatch (2000), P4-6 rat, 22 C.

Genesis implementation by Kai Du <kai.du@ki.se>, MScell v9.5.

Revision by Robert Lindroos <robert.lindroos@ki.se>, q factor applied
to both m and h instead of h only.

NEURON implementation by Alexander Kozlov <akozlov@csc.kth.se>.

ENDCOMMENT
