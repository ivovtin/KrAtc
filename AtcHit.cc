#include <iostream>
#include <cmath>
#include <string>

#include "AtcHit.hh"

#include "AtcVFitter.h"
#include "AtcPar.hh"

using namespace std;

int AtcHit::digitMode = READ5;
float AtcHit::ticks[5] = {0.,1.,2.,3.,4.};

int AtcHit::badEvent = AtcHit::BadAmplitude|AtcHit::Underflow;

const char * const AtcHit::typeLabels[nTypes] = {                    //типы сигналов
"BadAmplitude",	"Underflow",  	"Overflow",
"SomeAltPed",	"SomeAgtPed", 	"Pedestal",
"A0ltPed", 		"A0gtPed",		"NoPeak",
"OnePeak",		"TwoPeaks",		"Saw",
"ShiftedLeft",	"Tail",			"Rise",
"Spike",		"TypeDet",		"Estimated",
"TriedToFit",	"Fitted" };

AtcVFitter* AtcHit::leamaxFitter = 0;
AtcVFitter* AtcHit::minuitFitter = 0;
AtcVFitter* AtcHit::fitter = 0;

float AtcHit::maxAmpDevSig=3.;

float AtcHit::spikeMinDev=3., AtcHit::spikeMinRatio=2.5;

int AtcHit::verbose = 0;

AtcHit::AtcHit(int c) : cnt(c),
 A0(0.),  A(0.),  T(0.),  Tau(1.*tickSize),  P(5.), ChiSqr(0.),
dA0(0.), dA(0.), dT(0.), dTau(0.),          dP(0.), Tfinal(0.),
Nphe(0.), type(0), maxpos(0), peakmask(0), readpt(0), nread(0)
{
	a[0]=-1; a[1]=-1; a[2]=-1; a[3]=-1; a[4]=-1;
	da[0]=0.; da[1]=0.; da[2]=0.; da[3]=0.; da[4]=0.;
	err=1.;
}

AtcHit::AtcHit(int c,short *amp) : cnt(c),
 A0(0.),  A(0.),  T(0.),  Tau(1.*tickSize),  P(5.), ChiSqr(0.),
dA0(0.), dA(0.), dT(0.), dTau(0.),          dP(0.), Tfinal(0.),
Nphe(0.), type(0), maxpos(0), peakmask(0), readpt(0), nread(0)
{
	a[0]=amp[0]; a[1]=amp[1]; a[2]=amp[2]; a[3]=amp[3]; a[4]=amp[4];
	da[0]=0.; da[1]=0.; da[2]=0.; da[3]=0.; da[4]=0.;
	err=1.;
}

void AtcHit::setDigitMode(int mode)
{
	if(mode==READ4p1) {
		ticks[0]=0.;
        ticks[1]=3.*tickSize;
        ticks[2]=4.*tickSize;
        ticks[3]=5.*tickSize;
        ticks[4]=6.*tickSize;
	}
	else if(mode==READ5) {
        for(int i=0; i<5; i++)	ticks[i]=i*tickSize;
	}
	else {
		cerr<<"AtcHit::setDigitMode(): Unknown digitizing mode "<<mode<<endl;
		return;
	}
	digitMode=mode;
}

int AtcHit::init()
{
	if( cntPar[cnt].getfA0()>0 )	A0=cntPar[cnt].getfA0();
	if( cntPar[cnt].getError()>0 ) err=cntPar[cnt].getError();
	if( cntPar[cnt].getfTau()>0 )  Tau=cntPar[cnt].getfTau();
	if( cntPar[cnt].getfP()>0 )      P=cntPar[cnt].getfP();

	if( A0==0.0 ) {
	// Pedestal is undefined
		if( a[0]<0 ) return -1;
		A0=a[0];
		for(int i=0; i<5; i++) {
			if( a[i] >= 0 ) {
				readpt|=1<<i;
				nread++;
				if( da[i]==0.0 ) da[i] = err;
			}
		}
		type |= Pedestal;
	} else {
		float pedUpLim=A0+maxAmpDevSig*err;
		for(int i=0; i<5; i++) {
			if(  a[i] < 0 ) { //conditional reading
				a[i] = (int)A0;
				if( da[i]==0.0 ) da[i] = 2.;
			} else {
				readpt|=1<<i;
				nread++;
				if( da[i]==0.0 ) da[i] = err;
			}
		}
	}
	return 0;
}

int AtcHit::check4p1()
{
	int i;

	float pedLoLim=A0-maxAmpDevSig*err, pedUpLim=A0+maxAmpDevSig*err;

	badAmp=0;
	type=0;
	peakmask=0;
	int nPeaks=0;

	for(i=0; i<5; i++)
	{
		short ac=a[i], ap=a[i==0?0:i-1], an=a[i==4?4:i+1];
		uint8_t pointmask = 1<<i;
		if( ac>pedUpLim ) {
			type |= SomeAgtPed;
			if( i==0 ) type |= A0gtPed;
			if( ac==1023 ) { type |= Overflow; badAmp|=pointmask; }
			else if( ac>1023 ) { type |= BadAmplitude; badAmp|=pointmask; }
		} else if( ac<pedLoLim ) {
			type |= SomeAltPed;
			if( i==0 ) type |= A0ltPed;
			if( ac==0 ) { type |= Underflow; badAmp|=pointmask; }
			else if( ac<0 ) { type |= BadAmplitude; badAmp|=pointmask; }
		}
		if( (i==0 || ac>ap) && (i==4 || ac>=an) && ac>pedUpLim ) {
			peakmask|=pointmask;
			if( !nPeaks || ac>a[maxpos] ) maxpos=i;
			nPeaks++;
		}
		if( badAmp&pointmask && ac<=0 ) da[i]=A0;
	}

	if( !(type&SomeAgtPed || type&SomeAltPed) ) type |= Pedestal;


	if( nPeaks ) {
		switch( nPeaks ) {
			case 1: type |= OnePeak; break;
			case 2: type |= TwoPeaks; break;
			case 3: type |= Saw; break;
		}
		if( maxpos==0 ) {
			if( a[0]==a[1] ) maxpos++;
			else type |= Tail;
		}
		if( maxpos==1 ) {
			if( a[1]==a[2] ) maxpos++;
			else type |= ShiftedLeft;
		}
		if( maxpos==4 ) type |= Rise;
	}
	else
		type |= NoPeak;

	type |= TypeDet;

	return type;
}

int AtcHit::check5()
{
	int i;

    //limits of pedestal values
	float pedLoLim=A0-maxAmpDevSig*err, pedUpLim=A0+maxAmpDevSig*err;

	badAmp=0;
	type=0;
	peakmask=0;
	int nPeaks=0;

	for(i=0; i<5; i++)
	{
		short ac=a[i], ap=a[i==0?0:i-1], an=a[i==4?4:i+1];
		uint8_t pointmask = 1<<i;
		if( ac>pedUpLim ) {
			type |= SomeAgtPed;
			if( ac==1023 ) { type |= Overflow; badAmp|=pointmask; }
			else if( ac>1023 ) { type |= BadAmplitude; badAmp|=pointmask; }
		} else if( ac<pedLoLim ) {
			type |= SomeAltPed;
			if( ac==0 ) { type |= Underflow; badAmp|=pointmask; }
			else if( ac<0 ) { type |= BadAmplitude; badAmp|=pointmask; }
		}
		if( (i==0 || ac>ap) && (i==4 || ac>=an) && ac>pedUpLim ) {
			peakmask|=pointmask;
			if( !nPeaks || ac>a[maxpos] ) maxpos=i;
			nPeaks++;
		}
		if( badAmp&pointmask && ac<=0 ) da[i]=A0;
	}

	if( !(type&SomeAgtPed || type&SomeAltPed) ) type |= Pedestal;

	if( nPeaks ) {
		switch( nPeaks ) {
			case 1: type |= OnePeak; break;
			case 2: type |= TwoPeaks; break;
			case 3: type |= Saw; break;
		}
		if( maxpos==0 ) {
			if( a[0]==a[1] ) maxpos++;
			else type |= Tail;
		}
		if( maxpos==4 ) type |= Rise;

		//check the hit for spike
		short submax=0, submin=1024;

		for(int i=0; i<5; i++) {
			if( i!=maxpos ) {
				if( submax<a[i] ) submax=a[i];
				if( submin>a[i] ) submin=a[i];
			}
		}
		if( a[maxpos] > submax+spikeMinDev*err &&
			a[maxpos]-(submax+submin)/2 > spikeMinRatio*(submax-submin) )
			type|=Spike;
	} else
        type |= NoPeak;

	type |= TypeDet;

	return type;
}

int AtcHit::estimate()
{
	if( !type || type&NoPeak ) return 1;

	type|=Estimated;

	int i=maxpos;

	A=a[i]-A0;
	T=ticks[i];

	if( nread==1 ) return 1;
	if( type&(Rise|Tail|ShiftedLeft) ) return 0;

	float C=2*a[i]-a[i-1]-a[i+1], D=a[i+1]-a[i-1];

	if( C!=0 ) {
		A+=D*D/8/C;
		T+=tickSize*D/2/C;
	}

	return 0;
}

int AtcHit::rec(bool doNotFit)
{
	if( verbose ) cout<<"HIT RECONSTRUCTION IN COUNTER "<<cnt<<endl;

	if( init()!=0 ) {
		if( verbose ) cout<<"UNDEFINED PEDESTAL IN CONDITIONAL READING"<<endl;
		return -1;
	}

	if( type&Pedestal ) return 0;

	if(digitMode==READ4p1)
		check4p1();
	else
		check5();

	if( type&badEvent ) {
		if( verbose ) cout<<"BAD HIT DETERMINED"<<endl;
		return 1;
	}

	if( type&Pedestal ) {
		A=0;
		A0=a[0];
		return 0;
	}

	if( estimate()!=0 || doNotFit ) {
		if( cntPar[cnt].getA1pe() > 0 ) {
			Nphe=A/cntPar[cnt].getA1pe();               //вычисление числа ф.э. с использованием 1ф.э. амплитуды
			dNphe=0.1*Nphe; //estimate Nphe error as 10%
		} else { //set roughly Nphe even if A1pe not defined
			Nphe=A/30.;
			dNphe=5./30.*A; //error comes from +-5 counts possible A1pe variation
		}
		return 0;
	}

	type|=TriedToFit;

	if( verbose ) {
		cout.precision(4);
		cout<<"PARAMETERS ESTIMATION: "
			<<"A0="<<A0<<" A="<<A<<" T="<<T<<" Tau="<<Tau<<" P="<<P<<endl;
		cout.precision(0);
	}

	if( !fitter ) return -1; //no fitter declared though used

	int fitres = fitter->fit(this);

	if( fitres==0 ) type|=Fitted;

	if( cntPar[cnt].getA1pe() > 0 ) {
		Nphe=A/cntPar[cnt].getA1pe();
		dNphe=dA/cntPar[cnt].getA1pe();
	} else { //set roughly Nphe even if A1pe not defined
		Nphe=A/30.;
		dNphe=5./30.*Nphe; //error comes from +-5 counts possible A1pe variation
	}
	if( verbose ) {
		if( fitres!=0 ) {
			cout<<"FIT FAILED WITH RETURN STATUS "<<fitres<<". FINAL PARAMETER VALUES:\n";
			cout.precision(4);
			for(int i=0; i<fitter->getNpar(); i++)
				cout<<fitter->getParName(i)<<"="<<fitter->getParameter(i)<<" ";
			cout<<"ChiSqr/Ndf="<<fitter->getChiSqr()<<"/"<<fitter->getNdf()<<endl;
			cout.precision(0);
		}
	}
	return 0;
}

void AtcHit::print(ostream& out) const
{
	out<<"Counter "<<cnt<<"\nAmplitudes: ";
	for(int i=0; i<5;i++)
        out<<a[i]<<" ";
	out<<"("<<nread<<" read)"<<endl;

	if( type ) {
		out<<"Type: ";
	    for(int i=0; i<nTypes; i++)
			if( type&1<<i ) out<<typeLabels[i]<<" ";
	    out<<"  Peaks: ";
	    for(int i=0; i<5; i++)
    	   	if( peakmask&1<<i ) out<<i;
		out<<"  for fit "<<(int)maxpos<<endl;
		out.precision(4);
		out<<"A0="<<A0<<" A="<<A<<" T="<<T;
		if( Tfinal!=0. ) out<<" Tcorr="<<Tfinal;
		if( type&Fitted ) {
			out<<" Tau="<<Tau<<" P="<<P<<" ChiSqr/Ndf="<<ChiSqr<<"/"<<Ndf;
			if( fitter==leamaxFitter )
                out<<" (LEAMAX)";
			else if( fitter==minuitFitter )
                out<<" (Minuit)";
		}
        if( cntPar[cnt].getA1pe()>0 )
			out<<" Npe="<<Nphe;
        else
			out<<" Npe=undefined";
		out<<endl;
		out.precision(0);
	}
}

#ifndef ATC_NOROOT
# include "A6PulseShape.h"
# include "TMath.h"
# include "TPad.h"
# include "TStyle.h"
# include "TGraphErrors.h"
# include "TPolyMarker.h"
# include "TH1.h"
# include "TF1.h"
# include "TLegend.h"
# include "TList.h"

/*
Options:
'v' - print hit to stdout
's' - scale vertical axis proportional to hit amplitude
'c' - calibration run assumed
't' - draw pulse shape even if it is not fitted
'f' - do fit with another fitter and plot results
*/
int AtcHit::draw(TPad *pad,const char* opt) const
{
	static char title[100];
	static TGraphErrors gpread, gpunread;
	static TF1 *fFitPulse = 0, *fFitPulse2 = 0;
	static TH1S *hf=0;
	static TLegend *leg=0;
    static const Int_t kNotDraw=1<<9;

	if( !type ) return -1;

	string sopt(opt);

	if( sopt.find('v')!=string::npos ) print(cout);

//	cout.precision(3);
//	cout<<"Counter "<<cnt<<": "<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<" "<<a[4]<<endl;

	if( !hf ) {
		hf=new TH1S("hframe","",10,ticks[0]-1,ticks[4]+1);
		hf->SetMinimum(0);
		hf->SetStats(0);

		fFitPulse = new TF1("fFitPulse",A6PulseShape,-10,10,5);
		fFitPulse->SetLineColor(kBlue);
		fFitPulse->SetBit(kNotDraw,kTRUE);
		for(int i=0; i<fitter->getNpar(); i++)
			fFitPulse->SetParName(i,fitter->getParName(i));

		fFitPulse2 = new TF1("fFitPulse2",A6PulseShape,-10,10,5);
		fFitPulse2->SetLineColor(kGreen+2);
		fFitPulse2->SetBit(kNotDraw,kTRUE);
		for(int i=0; i<fitter->getNpar(); i++)
			fFitPulse2->SetParName(i,fitter->getParName(i));

		leg = new TLegend(0.15,0.73,0.35,0.88);
        leg->SetBorderSize(1);
		leg->SetFillStyle(0);

		//gpread.GetListOfFunctions()->Add(fFitPulse);
		//gpread.GetListOfFunctions()->Add(fFitPulse2);
		gpread.SetMarkerStyle(20);
		gpunread.SetMarkerStyle(24);

		gStyle->SetOptFit(1111);
		gStyle->SetStatStyle(0);
	}

    gpread.Set(nread);
	gpunread.Set(5-nread);
    int ipr=0, ipu=0;
	for(int i=0; i<5; i++) {
		if( isPointRead(i) ) {
			gpread.SetPoint(ipr,ticks[i],a[i]);
			gpread.SetPointError(ipr,0,da[i]);
			ipr++;
		} else {
			gpunread.SetPoint(ipu,ticks[i],a[i]);
			gpunread.SetPointError(ipu,0,da[i]);
			ipu++;
		}
	}

	if( sopt.find('c')!=string::npos )
		sprintf(title,"Counter %d  A=%.1f;Time, ticks;Amplitude, counts",cnt,A);
	else
		sprintf(title,"Counter %d  Npe=%.1f;Time, ticks;Amplitude, counts",cnt,Nphe);

	hf->SetTitle(title);

	pad->Clear();
    leg->Clear();
	fFitPulse->SetBit(kNotDraw,kTRUE);
	fFitPulse2->SetBit(kNotDraw,kTRUE);

	hf->Draw("0");
	gpread.Draw("p1");
	leg->AddEntry(&gpread,"Read points","p");
	if( nread<5 ) {
		gpunread.Draw("p1");
		leg->AddEntry(&gpunread,"Unread points","p");
	}

	if( sopt.find('s')!=string::npos )
		hf->SetMaximum(1.1*TMath::MaxElement(5,a));
    else
		hf->SetMaximum(1200);

	if( sopt.find('t')!=string::npos || (type&Fitted) ) {
		fFitPulse->ResetBit(kNotDraw);
		fFitPulse->SetParameters(A0,A,T,Tau,P);
		fFitPulse->SetParError(0,0);
		fFitPulse->SetParError(1,dA);
		fFitPulse->SetParError(2,dT);
		fFitPulse->SetParError(3,dTau);
		fFitPulse->SetParError(4,dP);
		fFitPulse->SetChisquare(ChiSqr);
		fFitPulse->SetNDF(Ndf);
//		cout<<"A0="<<A0<<" A="<<A<<" T="<<T<<" Tau="<<Tau<<" P="<<P<<" Chi2/ndf="<<ChiSqr<<'/'
//			<<ndf<<" Npe="<<Nphe<<endl;
		if( fitter==leamaxFitter )
	        leg->AddEntry(fFitPulse,"LEAMAX fit","l");
		else if( fitter==minuitFitter )
			leg->AddEntry(fFitPulse,"Minuit fit","l");
	}

	if( sopt.find('f')!=string::npos && (fitter==leamaxFitter || fitter==minuitFitter) && (type&TriedToFit) ) {
        AtcVFitter* oldFitter=fitter;
		if( fitter==leamaxFitter )
			fitter=minuitFitter;
		else
			fitter=leamaxFitter;

        AtcHit ahit(cnt);

		for(int i=0; i<5; i++) {
			if( isPointRead(i) ) ahit.a[i] = a[i];
		}

		ahit.rec();

		if( sopt.find('t') || ahit.checkType(Fitted) ) {
			fFitPulse2->ResetBit(kNotDraw);
			fFitPulse2->SetParameters(ahit.A0,ahit.A,ahit.T,ahit.Tau,ahit.P);
			fFitPulse2->SetParError(0,0);
			fFitPulse2->SetParError(1,ahit.dA);
			fFitPulse2->SetParError(2,ahit.dT);
			fFitPulse2->SetParError(3,ahit.dTau);
			fFitPulse2->SetParError(4,ahit.dP);
			fFitPulse2->SetChisquare(ahit.ChiSqr);
			fFitPulse2->SetNDF(ahit.Ndf);
			cout<<">>>Comparison with "<<(fitter==leamaxFitter?"LEAMAX":"Minuit")<<" fit:"
				<<" A1/A2="<<A/ahit.A<<" T1-T2="<<T-ahit.T<<" Tau1/Tau2="<<T/ahit.T<<" ChiSqr2/Ndf="<<ahit.ChiSqr<<"/"<<ahit.Ndf<<endl;
			if( fitter==leamaxFitter )
		        leg->AddEntry(fFitPulse2,"LEAMAX fit","l");
			else
				leg->AddEntry(fFitPulse2,"Minuit fit","l");
			leg->Draw();
		}

		fitter=oldFitter;
	}

	leg->Draw();
	pad->Modified();
	pad->Update();

	return 0;
}
#else //provide dummy draw() method if visualization disabled
int AtcHit::draw(TPad *pad,const char* opt) const
{
	cerr<<"ATC visualization disabled in library."<<endl;
	return 0;
}
#endif

