#ifndef AtcLeamaxFitter_hh
# define AtcLeamaxFitter_hh 1

# include "AtcVFitter.h"
# include "A6PulseShape.h"

class AtcHit;

class AtcLeamaxFitter : public AtcVFitter {
private:
	double parl[5], paru[5]; //lower and upper bounds of parameters
	bool fixed[5]; //mask of fixed parameters
	int verbose; //verbosity level

public:
	AtcLeamaxFitter(int verb=0);
	virtual ~AtcLeamaxFitter() {}

	void fix(int i) { if(i<5&&i>=0) fixed[i]=true; }
	void release(int i) { if(i<5&&i>=0) fixed[i]=false; }
	void releaseAll() { for(int i=0; i<5; i++) fixed[i]=false; }
	bool isFixed(int i) const { return fixed[i]?true:false; }

	int fit(AtcHit* aHit);

	double eval(double x) const { return A6PulseShape(x,par); }
	double eval(double x,double *p) const { return A6PulseShape(x,p); }

	void setVerboseLevel(int v) { verbose=v; }

private:
	static void sub(int *k,double x[],int *n,double a[],double *f,double df[],int *mode,int *nerror);
};

#endif
