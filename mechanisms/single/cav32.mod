TITLE T-type calcium channel (Cav3.2)

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
    SUFFIX cav32
    USEION cal READ cali, calo WRITE ical VALENCE 2
    RANGE pbar, ical, a, perm, I
}

PARAMETER {
    pbar = 6.7e-6   (cm/s)
    mvhalf = -61.5  (mV)
    mslope = -8.0   (mV)
    hvhalf = -73.7  (mV)
    hslope = 9.1   (mV) :9.1 
    a      = 0.9
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
    mtau  (ms)
    htau  (ms)
    htau2 (ms)
    htot  (ms)
    perm
    I
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    perm = pbar*m*m*m*h
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
    minf  = 1/(1+exp((v-mvhalf)/mslope))
    hinf  = 1/(1+exp((v-hvhalf)/hslope))    
    mtau  = 6.0/(1+exp((v+66.0)/15.0  ))+0.6
    htau  = 4.3/(1+exp(0.06*(v)))+8
    htau2 = 95*exp(-(v+58.0)/25.0)+20
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

Original data by Iftinca (2006) , rat, 37 C.

Genesis implementation by Kai Du (21 C) <kaidu828@gmail.com> m^2*h.

NEURON implementation by Alexander Kozlov <akozlov@nada.kth.se>, smooth
fit of mtau and htau.

Revised NEURON model by Robert Lindroos
-> 37 C (kinetics and infinity parameters)
-> Half activation and slope factors were found that approximated the activation curve 
   from the paper when cubed (m3). m3 were used in order to remove the large current 
   obtained when stepping back to holding potential (-110 mV), at the end of the trace. 
-> slow and fast inactivation was approximated from the paper and combined as
   tauh = 0.9 * fast + 0.1 * slow. Activation kinetics was also extracted from paper.

ENDCOMMENT
