#ifndef ATC_NOROOT

#include "AtcMinuitFitter.h"
#include "AtcHit.hh"
#include <iostream>
using std::cout;
using std::endl;

static const char* names[5]={ "A0", "A", "T", "Tau", "P" };

static AtcHit* Ihit = 0;
static void *libmin;

AtcMinuitFitter::AtcMinuitFitter(int verb) : AtcVFitter()
{
	Npar=5;
	ParNames = names;
	par=new double[Npar];
	dpar=new double[Npar];
	verbose=verb;

	atcMinuit = new TMinuit(Npar);
	atcMinuit->SetFCN(fcn);
	atcMinuit->SetPrintLevel(verbose-1); //set print out level for Minuit
	if( verbose==0 ) atcMinuit->Command("SET NOW"); //shut up warnings
	atcMinuit->Command("SET STRATEGY 0"); //minimum FCN calls
	atcMinuit->Command("SET GRA 1"); //FCN calculates gradient
	atcMinuit->SetErrorDef(1.0); //using chi^2 method
	atcMinuit->SetMaxIterations(200); //set max number of iteration
//	atcMinuit->Command("SET EPS 1.E-9"); // set machine accuracy

	atcMinuit->DefineParameter(ATCPAR_A0,ParNames[ATCPAR_A0],0.,1.,0.,0.);
	atcMinuit->DefineParameter(ATCPAR_A,ParNames[ATCPAR_A],0.,10.,0.,0.);
	atcMinuit->DefineParameter(ATCPAR_T,ParNames[ATCPAR_T],0.,AtcHit::tickSize,0.,0.);
	atcMinuit->DefineParameter(ATCPAR_TAU,ParNames[ATCPAR_TAU],1.*AtcHit::tickSize,AtcHit::tickSize,0.5*AtcHit::tickSize,2*AtcHit::tickSize);
	atcMinuit->DefineParameter(ATCPAR_P,ParNames[ATCPAR_P],2.,.1,0.1,10.);
	for(int i=0; i<5; i++) fixed[i]=false;
}
AtcMinuitFitter::~AtcMinuitFitter() {
	delete [] par;
	delete [] dpar;
	delete atcMinuit;
}
int AtcMinuitFitter::setPar(int idPar,float value)
{
	Double_t arglist[2]={idPar+1,value};
	Int_t ierflg;

	atcMinuit->mnexcm("SET PAR",arglist,2,ierflg);
	return ierflg;
}
void AtcMinuitFitter::fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
	Double_t chisq=0, a, delta, der[5], *df=0;
	Int_t i, k;
	if( iflag==2 ) {
		df=der;
		gin[0]=gin[1]=gin[2]=gin[3]=gin[4]=0.0;
	}
	for(i=0; i<5; i++) {
		a=A6PulseShape(AtcHit::ticks[i],par,df);
		delta = (Ihit->getamp(i)-a)/Ihit->geterra(i);
		if( iflag==2 )
			for(k=0; k<5; k++) gin[k]+=-2*der[k]*delta/Ihit->geterra(i);
		chisq += delta*delta;
	}
	f = chisq;
}
double AtcMinuitFitter::evalFCN(AtcHit* aHit,double *p)
{
	int n=5;
	double chi2;
	Ihit=aHit;
	fcn(n,0,chi2,p,0);
	return chi2;
}

int AtcMinuitFitter::fit(AtcHit* aHit)
{
	Ihit=aHit;

	setPar(ATCPAR_A0, aHit->getA0());
	setPar(ATCPAR_A,  aHit->getA());
	setPar(ATCPAR_T,  aHit->getT());
	setPar(ATCPAR_TAU,aHit->getTau());
	setPar(ATCPAR_P,  aHit->getP());
	bool releaseTau=false;
	if( aHit->getNread()<=2 && !fixed[ATCPAR_TAU] ) { //scarce data in conditional reading
		atcMinuit->FixParameter(ATCPAR_TAU); //fix tau for this hit only
		releaseTau=true;
	}

	int rc = atcMinuit->Migrad();

	atcMinuit->GetParameter(ATCPAR_A0, par[0],dpar[0]);
	atcMinuit->GetParameter(ATCPAR_A,  par[1],dpar[1]);
	atcMinuit->GetParameter(ATCPAR_T,  par[2],dpar[2]);
	atcMinuit->GetParameter(ATCPAR_TAU,par[3],dpar[3]);
	atcMinuit->GetParameter(ATCPAR_P,  par[4],dpar[4]);
	chi2=evalFCN(aHit,par);
	ndf=5-atcMinuit->GetNumFreePars();

	if( rc==0 ) {
		if( par[1]>=0.0 && par[1]<1e4 &&
			par[2]>-10*AtcHit::tickSize && par[2]<10*AtcHit::tickSize ) {
			aHit->setA0(par[0],dpar[0]);
			aHit->setA(par[1],dpar[1]);
			aHit->setT(par[2],dpar[2]);
			aHit->setTau(par[3],dpar[3]);
			aHit->setP(par[4],dpar[4]);
			aHit->setChiSqr(chi2);
			aHit->setNdf(ndf);
		} else { //wrong amplitude or time
			rc=10;
		}
	}

	if( releaseTau )
		atcMinuit->Release(ATCPAR_TAU); //release tau

	return rc;
}

#endif //ATC_NOROOT

