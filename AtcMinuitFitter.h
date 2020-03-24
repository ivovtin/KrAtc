#ifndef AtcMinuitFitter_hh
# define AtcMinuitFitter_hh 1

# include "AtcVFitter.h"
# include "Rtypes.h"
# include "TMinuit.h"
# include "A6PulseShape.h"

class AtcHit;

class AtcMinuitFitter : public AtcVFitter {
private:
	bool fixed[5];
	int verbose; //verbosity level
    TMinuit *atcMinuit;

public:
	AtcMinuitFitter(int verb=0);
	virtual ~AtcMinuitFitter();

	void fix(int i) { if(i<5&&i>=0) { fixed[i]=true; atcMinuit->FixParameter(i); } }
	void release(int i) { if(i<5&&i>=0) { fixed[i]=false; atcMinuit->Release(i); } }
	void releaseAll() { for(int i=0; i<5; i++) fixed[i]=false; atcMinuit->mnfree(0); } //release all currently fixed parameters

	int fit(AtcHit* aHit);

	double eval(double x) const { return A6PulseShape(x,par); }
	double eval(double x,double *p) const { return A6PulseShape(x,p); }
	static double evalFCN(AtcHit* aHit,double *p);

	void setVerboseLevel(int v) { verbose=v; }

private:
	int setPar(int idPar,float value); //set initial values of Minuit parameters

	static void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
};

#endif
