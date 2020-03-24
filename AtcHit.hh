#ifndef AtcHit_hh
# define AtcHit_hh

# include <stdint.h>
# include <iostream>
# include <map>

class AtcVFitter;
class TPad;

//***************************************************************
//* AtcHit stores event information of one counter              *
//* as well as reconstructed amplitude and time, event type, etc*
//***************************************************************
class AtcHit {
    friend class AtcRec;
public:
	enum EventType {  //to be determined before fit
		BadAmplitude  =1,       //some amplitude absent or >1023
		Underflow     =2,       //some amplitude is 0
		Overflow      =1<<2,    //some amplitude is 1023
		SomeAltPed    =1<<3,    //some amplitude is less than pedestal
		SomeAgtPed    =1<<4,    //some amplitude is grater than pedestal
		Pedestal      =1<<5,    //all amplitudes are near pedestal
		A0ltPed       =1<<6,    //if first point less than pedestal
		A0gtPed       =1<<7,    //if first point grater than pedestal
		NoPeak        =1<<8,    //no peak found
		OnePeak       =1<<9,    //one peak found
		TwoPeaks      =1<<10,   //two peaks found
		Saw           =1<<11,   //three peaks found
		ShiftedLeft   =1<<12,   //maximum peak is in point 2, old digitizing only
		Tail          =1<<13,   //maximum peak is in point 1
		Rise          =1<<14,   //maximum peak is in point 5
		Spike         =1<<15,   //weird spike in a single point (observed in calibration with generator)
		TypeDet       =1<<16,   //type determining member function check() called
		Estimated     =1<<17,   //amplitude estimating function estimate() called
		TriedToFit    =1<<18,   //fit was attempted
		Fitted        =1<<19   //fit success returned by Minuit
	};

	static const int nTypes=20;

	static const char * const typeLabels[nTypes];

	static const float tickSize=1.; //distance between ticks
	static float ticks[5]; //A6 time ticks

	static float maxAmpDevSig; //acceptable amplitude deviation in units of sigmas

	//parameters for determining spikes:
	// spikeMinDev    - deviation of peak amplitude from substrate in number of errors
	// spikeMinRatio  - minimum ratio of peak amplitude to substrate span
	static float spikeMinDev, spikeMinRatio;

	static int badEvent; //Bit mask of event types not fitted (see EventType)

	enum { READ4p1, READ5 }; //digitizing modes

	static AtcVFitter *leamaxFitter, *minuitFitter; //pointers to hit fitting managers
	static AtcVFitter *fitter; //fitter instance used to fit pulse shape

	static int verbose; //verbosity flag
private:
	int cnt; //Counter number 1..160

	short a[5]; //A6 amplitudes

	float da[5]; //A6 amplitude errors, used in fitting
	float err; //amplitude error of the channel

	float A0, A, T, Tau, P; //Fit parameters
	float ChiSqr; //Fit quality - Chi square
	int Ndf; //number of freedom degrees
	float dA0, dA, dT, dTau, dP; //Fit parameter errors
	float Tfinal; //Corrected to deltaT time

	float Nphe, dNphe; //Number of photoelectrons and its statistical error

	int type;         //Bit mask of event type determined by check()
	uint8_t badAmp;   //Bit mask of bad amplitudes determined by check()
	uint8_t maxpos, peakmask; //peak position for fit and peak mask determined by check()
	uint8_t readpt; //read points mask (0 - if digit is absent, 1 - otherwise)
	int nread; //number of points read

	static int digitMode; //current digitizing mode
public:
	AtcHit(int i=0);
	AtcHit(int i,short *amp);
	~AtcHit() {}

protected:
	int init(); //preset pedestal value, errors, skipped amplitudes
	int check4p1(); //check event type, get peak positions, set errors (for 4+1 mode)
	int check5(); //check for 5-points mode
	int estimate(); //quadratic estimation of peak amplitude and time

public:
	static void setDigitMode(int mode);

	int rec(bool doNotFit=false); //hit reconstruction

	void print(std::ostream& out=std::cout) const; //print hit parameters
	int draw(TPad *pad,const char* opt="vt") const; //draw hit

	void setamp(short *amp) { a[0]=amp[0]; a[1]=amp[1]; a[2]=amp[2]; a[3]=amp[3]; a[4]=amp[4]; }
	void seterr(float e) { err=e; }

	int getCounter() const { return cnt; }
	int getType() 	const { return type; }
	bool checkType(int mask) const { return (type&mask)?true:false; }
	bool checkTypeStrict(int mask) const { return (type&mask)==mask?true:false; }
	bool notType(int mask) const { return (type&mask)==0?true:false; }
	int getMaxPos() const { return maxpos; }
	uint8_t getReadPtMask() const { return readpt; }
	bool isPointRead(int i) const { return readpt&(1<<i)?true:false; }
	int   getNread() const { return nread; }
	short getamp(int i) const { return a[i]; }
	float geterra(int i) const { return da[i]; }
	float geterr() const { return err; }
	float getA0() 	const { return A0; }
	float getdA0() 	const { return dA0; }
	float getA() 	const { return A; }
	float getdA() 	const { return dA; }
	float getT() 	const { return T; }
	float getdT() 	const { return dT; }
	float getTau() 	const { return Tau; }
	float getdTau() const { return dTau; }
	float getP() 	const { return P; }
	float getdP() 	const { return dP; }
	float getChiSqr() const { return ChiSqr; }
	int   getNdf()    const { return Ndf; }
	float getTfinal() const { return Tfinal; }
	float getNphe()   const { return Nphe; }
	float getdNphe()  const { return dNphe; }

	void setA0(float v,float dv) 	{ A0=v; dA0=dv; }
	void setA(float v,float dv) 	{ A=v; dA=dv; }
	void setT(float v,float dv) 	{ T=v; dT=dv; }
	void setTau(float v,float dv) 	{ Tau=v; dTau=dv; }
	void setP(float v,float dv)		{ P=v; dP=dv; }
	void setChiSqr(float chi2)  	{ ChiSqr=chi2; }
	void setNdf(int ndf)	  		{ Ndf=ndf; }
};

class AtcHits : public std::map<int,AtcHit>              //доступ к объектам AtcHit через контейнер AtcHits
{
public:
	typedef std::map<int,AtcHit> base_type;
public:
	short deltaT;       //digit from DeltaT card
	float deltaTfinal;  //normalized deltaT
	int ev_type;        //type of event (0-OK, 1-spoilt raw record, 2-DeltaT out of range)

	AtcHits() : base_type(), deltaT(0), deltaTfinal(0.0), ev_type(0) {}

	void clear() {
		base_type::clear();
		deltaT=0;
		deltaTfinal=0.0;
		ev_type=0;
	}
};

#endif
