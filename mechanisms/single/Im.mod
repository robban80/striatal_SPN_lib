COMMENT
Mechanism taken from Doron et al., 2017
https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=231427&file=/reproduction/Im.mod#tabs-2

Reference :     Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones

corrected rates using q10 = 2.3, target temperature 34, orginal 21

---------------------------------------------------------------

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

NEURON	{
	SUFFIX Im
	USEION k READ ek WRITE ik
	RANGE gbar, gIm, ik
    RANGE damod, maxMod, level, max2, lev2
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2) 
    damod = 0
    maxMod = 1
    level = 0
    max2 = 1
    lev2 = 0
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gIm	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gIm = gbar*m*modulation()
	ik = gIm*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)

	UNITSOFF
		mAlpha = 3.3e-3*exp(2.5*0.04*(v - -35))
		mBeta = 3.3e-3*exp(-2.5*0.04*(v - -35))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = (1/(mAlpha + mBeta))/qt
	UNITSON
}

FUNCTION modulation() {
    : returns modulation factor
    
    modulation = 1 + damod * ( (maxMod-1)*level + (max2-1)*lev2 ) 
    if (modulation < 0) {
        modulation = 0
    } 
}
