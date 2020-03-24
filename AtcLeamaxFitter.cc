#include <cmath>
#include <iostream>
#include <unistd.h>

#include "AtcLeamaxFitter.h"
#include "AtcHit.hh"
#include "dfunft.h"

static const char* names[5]={ "A0", "A", "T", "Tau", "P" };

AtcLeamaxFitter::AtcLeamaxFitter(int verb) : AtcVFitter()
{
	Npar=5;
	ParNames = names;
	par=new double[Npar];
	dpar=new double[Npar];
	verbose=verb;

	parl[ATCPAR_A0]=0.; paru[ATCPAR_A0]=1023.;
	parl[ATCPAR_A]=0.; paru[ATCPAR_A]=1e4;
	parl[ATCPAR_T]=-10*AtcHit::tickSize; paru[ATCPAR_T]=10*AtcHit::tickSize;
	parl[ATCPAR_TAU]=0.5*AtcHit::tickSize; paru[ATCPAR_TAU]=2*AtcHit::tickSize;
	parl[ATCPAR_P]=0.1; paru[ATCPAR_P]=10.;
	for(int i=0; i<5; i++) fixed[i]=false;
}

int AtcLeamaxFitter::fit(AtcHit* aHit)
{
	//input arguments
	static double eps=1e-8;
	static int k=1, m=5, n=5, mode=1, maxit=100, iprt;
	double x[5]={AtcHit::ticks[0],AtcHit::ticks[1],AtcHit::ticks[2],AtcHit::ticks[3],AtcHit::ticks[4]};
	double y[5]={aHit->getamp(0),aHit->getamp(1),aHit->getamp(2),aHit->getamp(3),aHit->getamp(4)};
	double sy[5]={aHit->geterra(0),aHit->geterra(1),aHit->geterra(2),aHit->geterra(3),aHit->geterra(4)};
	double al[5], au[5];

	//output arguments
	static int mfr, iafr[10], nerr;
	static double phi, dphi[5], cov[25], wrkspace[300];

	par[0]=aHit->getA0();
	par[1]=aHit->getA();
	par[2]=aHit->getT();
	par[3]=aHit->getTau();
	par[4]=aHit->getP();
	for(int i=0; i<5; i++) {
		if( fixed[i] ) al[i]=au[i]=par[i];
		else al[i]=parl[i], au[i]=paru[i];
	}
	if( aHit->getNread()<=2 ) //scarce data in conditional reading
		al[ATCPAR_TAU]=au[ATCPAR_TAU]=par[ATCPAR_TAU]; //fix tau for this hit only
	if( verbose>0 )
		iprt=verbose>1?-1:1;
	else
		iprt=0;

	dfunft_(sub,&k,&m,&n,&k,&n,x,y,sy,par,al,au,&mode,&eps,&maxit,&iprt,&mfr,iafr,&phi,dphi,cov,dpar,wrkspace,&nerr);
	chi2=2*phi;
	ndf=5-mfr;

	int rc=nerr;

	if( nerr==0 || nerr==4 ) { //nerr==4 often means that errors are not determined
		rc=0;
		aHit->setA0(par[0],dpar[0]);
		aHit->setA(par[1],dpar[1]);
		aHit->setT(par[2],dpar[2]);
		aHit->setTau(par[3],dpar[3]);
		aHit->setP(par[4],dpar[4]);
		aHit->setChiSqr(chi2);
		aHit->setNdf(ndf);
	}

	return rc;
}

void AtcLeamaxFitter::sub(int *k,double x[],int *n,double a[],double *f,double df[],int *mode,int *nerr)
{
	double *der=0;
	*nerr=0;

	if( *mode==1 ) der=df;

	*f=A6PulseShape(*x,a,der);
}

