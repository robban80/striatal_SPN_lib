TITLE T-type calcium channel (Cav3.3)

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

NEURON {
    SUFFIX cav33
    USEION cal READ cali, calo WRITE ical VALENCE 2
    RANGE pbar, ical, mvhalf, hvhalf, a, p, perm, I
}

PARAMETER {
    pbar = 6.7e-6 (cm/s)
    mvhalf = -73.5 (mV)     : 73.5 +/- 1.3
    mslope =  -4.4 (mV)
    hvhalf = -73.4 (mV)     : 73.4 +/- 2.5
    hslope =   5.6 (mV)
    a      = 0.9
    p      = 2

}

ASSIGNED { 
    v (mV)
    ical (mA/cm2)
    ecal (mV)
    celsius (degC)
    cali (mM)
    calo (mM)
    minf
    hinf
    mtau (ms)
    htau (ms)
    htau2 (ms)
    htot (ms)
    perm
    I
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    perm = pbar*(m^p)*h
    ical = ghk(v, cali, calo)*perm
    I    = ical
}

INITIAL {
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states { 
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htot
}

PROCEDURE rates(v (mV)) {
    minf = 1.0/(1+exp((v-mvhalf)/mslope))
    hinf = 1.0/(1+exp((v-hvhalf)/hslope))
    mtau = 6.0/(1+exp((v+73.0)/14.0))+0.8/(1+exp( (v-5.0)/3.0))+0.7 
    htau = 5.5/(1 + exp(0.06*(v)))+7
    htau2 = 90*exp(-(v+68.0)/40.0)+10
    htot  = a*htau + (1-a)*htau2
}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    if(z == 0) {
        z = z+1e-6
    }
    eco = co*(z)/(exp(z)-1)
    eci = ci*(-z)/(exp(-z)-1)
    ghk = (1e-3)*2*FARADAY*(eci-eco)
}

COMMENT 

Original data by Iftinca (2006), rat, 37 C

Genesis implementation (22 C) by Kai Du <kaidu828@gmail.com>, m*h, also corrected
for m^2*h.

NEURON implementation by Alexander Kozlov <akozlov@nada.kth.se>, smooth
fit of mtau and htau.

Revised model by Robert Lindroos 
-> 37 C (kinetics and infinity parameters; before room temp were used and adjusted by q-factor)
-> m2 used directly on reported values (m) without adjustment (see, cav32.mod).
   This gives good fit to experimental IV curve.
-> slow and fast inactivation combined as 0.9*fast + 0.1*slow

ENDCOMMENT
