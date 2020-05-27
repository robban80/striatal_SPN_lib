TITLE SK-type calcium activated K channel

UNITS {
    (molar) = (1/liter)
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (millimolar)
}

NEURON {
    SUFFIX sk
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gbar, ik, I
}

PARAMETER {
    gbar = 0.0 (mho/cm2)
    q = 1
}

ASSIGNED {
    v (mV)
    ik (mA/cm2)
    cai (mM) 
    ek (mV)
    oinf
    otau (ms)
    I (mA/cm2)
}

STATE { o }

BREAKPOINT {
    SOLVE state METHOD cnexp
    I  = gbar*o*(v-ek)
    ik = I 
}

DERIVATIVE state {
    rate(v, cai)
    o' = (oinf-o)/otau*q
}

INITIAL {
    rate(v, cai)
    o = oinf
}

PROCEDURE rate(v (mV), ca (mM)) {
    LOCAL a
    :ca = ca-(-0.0002) : Lindroos
    a = (ca/0.57e-3)^5.2
    oinf = a/(1+a)
    otau = 4.9
}

COMMENT

Original data by Hirschberg (1998) and Maylie (2004), Xenopus oocytes, room temp.

Genesis implementation by Evans (2013).

Revision by Robert Lindroos <robert.lindroos@ki.se>, Ca conc is shifted by -0.0002.
         done in order to shift the channel into a conducting state.
         similar effect can be obtained by a large increase of the maximal conductance 

NEURON implementation by Alexander Kozlov <akozlov@csc.kth.se>.


ENDCOMMENT
