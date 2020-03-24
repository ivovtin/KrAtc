#include <math.h>

double A6PulseShape(double x,double par[5],double *der=0)
{
	double &a0=par[0], &am=par[1], &tm=par[2], &tau=par[3], &p=par[4];

	double t0=tm-tau*log(p+1);
	if( x<=t0 ) {
		if( der!=0 ) {
			der[0]=1;
			der[1]=der[2]=der[3]=der[4]=0.0;
		}
		return a0;
	}

	double val, e, b, pw;
	e=exp((tm-x)/tau);
	b=(1+p-e)/p;
	pw=pow(b,p);
	val=a0 + am*e*pw;
	if( val>=1023 ) {
		if( der!=0 ) {
			der[0]=der[1]=der[2]=der[3]=der[4]=0.0;
		}
		return 1023;
	}

	if( der!=0 ) { //return to der partial derivatives
		double pwp=pw/b;
		der[0]=1;
		der[1]=e*pw;
		der[2]=am/tau*e*pwp*(b-e);
		der[3]=der[2]*(x-tm)/tau;
		der[4]=-am*e*pwp*(1-e)/p;
	}

	return val;
}

//for ROOT TF1 construction
double A6PulseShape(double *xx,double *par)
{
	double &x=*xx, &a0=par[0], &am=par[1], &tm=par[2], &tau=par[3], &p=par[4];

	if( x<=tm-tau*log(p+1) ) return a0;

	double val, e;
	e=exp((tm-x)/tau);
	val=a0 + am*e*pow((1+p-e)/p,p);;
	if( val>1023 ) return 1023;
	return val;
}


