#include <Riostream.h>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include"TROOT.h"
#include"TCanvas.h"
#include"TFile.h"
#include <iomanip>
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include <vector>

#include "KcSys/fortran.h"
#include "ReadNat/re_def.h"

extern "C" {

vector<float> pmt_vec;
vector<int> ch_vec;

int atc_pmt_signal_(float POI_atc, int J)
{
    //cout<<POI_atc<<"\t"<<J<<endl;
    pmt_vec.push_back(POI_atc);
    ch_vec.push_back(J);
}

int atc_pmt_signal_fill_(int cnt, float npe)
{
    TFile *fout=0;
    fout = new TFile(TString::Format("/spool/users/ovtin/test/atc_pmt_signal_sh_%d_%d_%f.root",cnt,kedrraw_.Header.Number,npe).Data(),"RECREATE");
    //cout<<"Event="<<kedrraw_.Header.Number<<endl;
    TProfile* pr=new TProfile("PMT:ch","PMT:ch",1024,0.0,1024.0,0,5);
    for( u_int k=0; k<pmt_vec.size(); k++ )
    {
        pr->Fill(ch_vec.at(k),pmt_vec.at(k));
	//cout<<pm_vec.at(k)<<"\t"<<ch_vec.at(k)<<endl;
    }
    fout->Write();
    fout->Close();
    pmt_vec.clear();
    ch_vec.clear();
}

}
