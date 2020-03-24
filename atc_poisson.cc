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

vector<float> poisson_vec;
vector<int> j_vec;

int atc_poisson_(float POI_atc, int J)
{
    //cout<<POI_atc<<"\t"<<J<<endl;
    poisson_vec.push_back(POI_atc);
    j_vec.push_back(J);
}

int atc_poisson_fill_(int cnt, float npe)
{
    TFile *fout=0;
    fout = new TFile(TString::Format("/spool/users/ovtin/test/atc_poisson_sh_%d_%d_%f.root",cnt,kedrraw_.Header.Number,npe).Data(),"RECREATE");
    //cout<<"Event="<<kedrraw_.Header.Number<<endl;
    TProfile* pr=new TProfile("POI:N","POI:N",50,-0.5,49.5,0,1);
    for( u_int k=0; k<poisson_vec.size(); k++ )
    {
        pr->Fill(j_vec.at(k),poisson_vec.at(k));
	//cout<<poisson_vec.at(k)<<"\t"<<j_vec.at(k)<<endl;
    }
    fout->Write();
    fout->Close();
    poisson_vec.clear();
    j_vec.clear();
}

}
