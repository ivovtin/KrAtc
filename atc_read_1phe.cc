//skript for get calibration spectrum from root-file
#include <Riostream.h>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include"TROOT.h"
#include"TCanvas.h"
#include"TFile.h"
#include <iomanip>
#include "TH2.h"

#include "KcSys/fortran.h"
#include "VDDCRec/kdcpar.h"
#include "VDDCRec/kglobparam.h"
#include "KDB/kdb.h"

extern "C" {

    TFile *f_1phe=0;
    TH1F *h=0;
    int ii=0;

    inline char* timestamp(time_t t)
    {
	static const char* time_fmt="%b %d %Y %H:%M:%S";
	static char strtime[50];

	struct tm* brtm=localtime(&t);
	strftime(strtime,50,time_fmt,brtm);

	return strtime;
    }

  int atc_read_1phe_(int ICNT, int k, float &bin)
{
    if( ii<1 ){
	cout<<"\t"<<"MC ATC system -> NSimRun="<<NSimRun<<endl;

	KDBconn* connection=kdb_open();
	if( !connection ) {
	    cerr<<"Can not establish connection to database"<<endl;
	    return -1;
	}

	time_t runTime; //start time of the current run
	//load calibration only if it is different from previous one
	time_t caltime;
	// Database tables info
	enum DBtableID     { Ped_ID=1302,  Gen_ID=1305,  LED_ID=1307,  DT_ID=1308 };
	enum DBtableLength { Ped_len=322,  Gen_len=1284, LED_len=1283, DT_len=2   };
	time_t last1peTime; //last begintime of 1pe calibration
	int ledDB[LED_len];
	static int ledCalN=0;


	runTime=kdb_run_get_begin_time(connection, NSimRun);
	cout<<"Begin time of "<<NSimRun<<": "<<timestamp(runTime)<<endl;

	if( runTime<=0 ) {
	    kdb_close(connection);
	}

	caltime=kdb_get_begin_time(connection,LED_ID,runTime);
	if( last1peTime==0 || last1peTime!=caltime ) {
	    if( kdb_read(connection, LED_ID, runTime, ledDB, LED_len) ) {
		last1peTime=caltime;
		ledCalN=ledDB[LED_len-1];
		cout<<"1ph.e. calibration "<<ledCalN<<" loaded ("<<timestamp(caltime)<<")"<<endl;
	    } else {
		cerr<<"Failed to get 1pe calibration"<<endl;
	    }
	} else {
	    cout<<"Use loaded 1pe calibration "<<ledCalN<<" ("<<timestamp(caltime)<<")"<<endl;
	}

	kdb_close(connection);

	f_1phe=TFile::Open(TString::Format("/spool/users/atc/ResLedCalib/calib_res_%02d.root", ledCalN).Data());
	if(f_1phe==0)
	{
	    printf("File not opened: Use file /spool/users/atc/ResLedCalib/calib_res_2197.root ... \n");
	    f_1phe=TFile::Open("/spool/users/atc/ResLedCalib/calib_res_2197.root");    //define calibration run
	} else {
	    printf("File opened successfully\n");
	}

    }
    ii++;

    if( k<2 ){
	TCanvas * c1 = (TCanvas*)f_1phe->Get(TString::Format("Ampl_%02d", ICNT).Data());
	h = (TH1F*)c1->GetPrimitive(TString::Format("Result_%02d", ICNT).Data());
	printf("ATC-> 1ph.e. amplitude is read for counter %d\r",ICNT); fflush(stdout);
    }

    int Nxbin=h->GetNbinsX();
    bin=h->GetBinContent(k);
 }
}
