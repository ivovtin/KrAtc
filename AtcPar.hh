#ifndef AtcPar_hh
#define AtcPar_hh

#include "AtcGlobals.h"

//***************************************************************
//* AtcCntPar stores counter's fixed parameter values           *
//***************************************************************

class AtcCntPar {
	friend class AtcRec;
private:
	int cnt; //number of counter
	bool off; //whether is excluded from reconstruction
	float aError; //amplitude measurement error
	float fA0, fTau, fP; //fixed values of parameters
	float rmsTau, rmsP; //known errors for fixed parameters
	float A1pe, rmsA1pe; //one photoelectron amplitude (from LED calibration)

	float tCutLeft, dtCutLeft;  // time and deltaT cut for left size events
	float tCutRight, dtCutRight;//    - | -      - | -     right   - | -

public:
	AtcCntPar(int c=0,float a0=-1.,float error=1.,float tau=1.,float p=5.,
		   float dtau=0.,float dp=0.) :
		cnt(c), aError(error), fA0(a0), fTau(tau), fP(p), rmsTau(dtau), rmsP(dp),
	    A1pe(0.), rmsA1pe(0.)
	{}

	~AtcCntPar() {}

	void setCounter(int c) { cnt=c; }
	void setError(float err) { aError=err; }; //set measurement error
	void setfA0(float val,float err=1.) { fA0=val, aError=err; } //set initial value and error of pedestal
	void setA1pe(float val,float err=10.) { A1pe=val, rmsA1pe=err; } //set one photoelectron amplitude (from LED calibration)
	void setfTau(float val,float err=0.) { fTau=val, rmsTau=err; } //set initial value and error  of time constant
	void setfP(float val,float err=0.) { fP=val, rmsP=err; } //set initial value and error of power parameter
	void setTdTcut(float tl, float dtl, float tr, float dtr)
	     { tCutLeft=tl; dtCutLeft=dtl; tCutRight=tr; dtCutRight=dtr; }

	bool  isOffline() const     { return off; }
	float getError() const      { return aError; }
	float getfA0() const        { return fA0; }
	float getfTau() const       { return fTau; }
	float getfP() const 	    { return fP; }
	float getA1pe() const 	    { return A1pe; }
	float getA1peErr() const    { return rmsA1pe; }
	float getTcutLeft() const   { return tCutLeft; }
	float getTcutRight() const  { return tCutRight; }
	float getDTcutLeft() const  { return dtCutLeft; }
	float getDTcutRight() const { return dtCutRight; }
};

struct AtcPar {
	AtcCntPar cpar[ATC_NCNT+1];
	float minDeltaT, maxDeltaT;
};

extern AtcPar atcPar;

static AtcCntPar* cntPar=atcPar.cpar;

#endif
