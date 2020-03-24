#ifndef AtcVFitter_hh
# define AtcVFitter_hh 1

# include <cmath>

class AtcHit;

class AtcVFitter {
protected:
	int Npar;
	const char** ParNames; //Parameter names
	double *par, *dpar;
	double chi2;
	int ndf;

public:

	AtcVFitter() : Npar(0), ParNames(0), par(0), dpar(0), chi2(0), ndf(0) {}
	virtual ~AtcVFitter() {}

	virtual void fix(int i) = 0;
	virtual void release(int i) = 0;
	virtual void releaseAll() = 0;

	virtual int fit(AtcHit* aHit) = 0;

	virtual double eval(double x) const = 0;
	virtual double eval(double x,double *p) const = 0;

	int getNpar() const { return Npar; }
	const char* getParName(int i) const { return ParNames[i]; }
	double getParameter(int i) const { return (i<Npar&&i>=0?par[i]:NAN); }
	double getParameterError(int i) const { return (i<Npar&&i>=0?dpar[i]:NAN); }
	double getChiSqr() const { return chi2; }
	int getNdf() const { return ndf; }

	void initParameter(int i,double val) { if(i<Npar&&i>=0) par[i]=val; }
};

#endif
