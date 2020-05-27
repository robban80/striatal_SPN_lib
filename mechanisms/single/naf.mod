TITLE Fast transient sodium current

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
    SUFFIX naf
    USEION na READ ena WRITE ina
    RANGE gbar, gna, ina, mVhalf, hVhalf, mSlope, hSlope, taum, tauh, taun, mtau,htau
    RANGE damod, maxMod, level, max2, lev2
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 1.8
    mVhalf     = -25.0 (mV)
    hVhalf     = -62.0 (mV)
    mSlope     =  -9.2 (mV)
    hSlope     =   6.0 (mV)
    taum       =   0.09 (ms)
    taun       =   0.34 (ms)
    tauh       =   0.34 (ms)
    damod = 0
    maxMod = 1
    level = 0
    max2 = 1
    lev2 = 0
} 

ASSIGNED {
    v (mV)
    ena (mV)
    ina (mA/cm2)
    gna (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gbar*m*m*m*h*modulation()
    ina = gna*(v-ena)
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
    UNITSOFF
    minf = 1 / (1 + exp( (v-mVhalf) / mSlope ) )
    hinf = 1 / (1 + exp( (v-hVhalf) / hSlope ) )
    
    mtau = 0.38 + 1/( 0.6*exp((v-(-58.0))/8.0) + 1.8*exp((v-(-58.0))/(-35.0))  )
    
    if (v < - 60) {
        htau = 3.4 + 0.015*v
    }else{
        htau = 0.56 + 1.1/(1+exp((v-(-48))/15.0)) + 1.2/(1+exp((v-(-48))/4.0))
    }
        
    :mtau = 0.13 +1/(0.6*exp((v-(-58))/taum)+1.8*exp((v-(-58))/(taun)))
    :htau = 0.14 +1.2/(1+exp((v-(-32))/tauh))
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

TODO:
update kintetics of inactivation gate below -60 mV, from linear to sigmoidal?
--------------------------------------------------------------------------------------

Original data by Ogata (1990), guinea pig, 22 C.

Genesis implementation by Kai Du <kai.du@ki.se>, MScell v9.5.

NEURON implementation by Alexander Kozlov <akozlov@csc.kth.se>, smooth
fit of mtau and htau. 

Updates by Robert Lindroos
channel kinetics updated so that natural instead of base 10 logarithm is used for 
values in Ogata 1990. This also reduces the AHP and spike frequency.
The revision of Ogata 1990, was done following a scanning of m and h gate parameters
that reduced the AHP magnitude.

Q factor of 1.8 used. Based on 
"The effect of temperature on Na currents in rat myelinated nerve fibres"
-> Q10 between 40 and 20 C 1.8-2.1
lowest value chousen since spiking is still intact and this make the cell fire slightly slower

ENDCOMMENT
