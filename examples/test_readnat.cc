using namespace std;
#include <unistd.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <math.h>
#include "AppFramework/AppFramework.hh"
#include "ReadNat/ReadNat.hh"
#include "ReadNat/read_nat_c.h"
#include "KaFramework/mknatfilelist.h"

#include "KrAtc/AtcRec.hh"
#include "KdATCCalib/AtcLed.hh"
#include "KrAtc/AtcHit.hh"
#include "KrAtc/AtcPar.hh"
#include "KrAtc/DCTrackExt.h"                    //����������� ����������� � ��������� � ����������� ����� (������ � �������)
#include "KrAtc/atc_geometry.h"            //��������� ���������

/*
#include "VDDCRec/kdcpar.h"
#include "VDDCRec/kvdpar.h"
#include "VDDCRec/ktracks.h"
#include "VDDCRec/VDDCRec.hh"
#include "VDDCRec/kdcswitches.h"
#include "VDDCRec/kdcvd.h"
#include "VDDCRec/khits.h"
#include "VDDCRec/ToFTrack.hh"
extern "C" {
#include "VDDCRec/ktracks.h"
#include "VDDCRec/mtofhits.h"
}
*/
//#include "KEmcRec/emc_system.h"
//#include "KEmcRec/emc_struct.h"
//#include "KrToF/tof_system.h"
//#include "KrAtc/atcrec.h"
//#include "KrVDDCMu/dcmu.h"
//#include "KrMu/mu_system.h"
//#include "KrMu/mu_event.h"
/*
extern "C" {
   void kedr_open_nat(char *fname,int *err);
   void kedr_read_nat(int *err);
   void kedr_close_nat(int *err);
}
*/
//  ������� ��������� � ������������ ������ AppFramework � ������� ��� ����-������� addModule
AppFramework::AppFramework(const char* sequence){
	addModule( new ReadNat     ("ReadNat", "Data reading"              ) );       //���������� ������ ReadNat
//	addModule( new VDDCRec     ("VDDCRec", "Tracks reconstruction"     ) );
	addModule( new AtcRec     ("AtcRec", "Atc reconstruction"     ) );
////	addModule( new EmcRec     ("EmcRec", "Emc reconstruction"     ) );
////	addModule( new RecWrapper  ("RecWrapper", "root reconstruction"    ) );
	//	reorder(sequence);             //������� ������������������ �������
}

int main(int argc, char* argv[])
{

    AppFramework* fw=new AppFramework("ReadNat->AtcRec");    //AppFramework- ����-�� ���������� ����� ��������
    fw->modify("ReadNat::input","/home/ovtin/simulation/simout/sim000001.dat");      //�������� ������
    //fw->modify("ReadNat::restoreReco",1);                    //������ �������������
    //������� root-���� ������������� ��� ������������ ������� - �� ��������
    //fw->modify("RecWrapper::browsable",0);
    //fw->modify("RecWrapper::output","./output.root");
    //fw->modify("RecWrapper::content","saveRawHeader:saveRawVD:saveRawDC:saveRawToF:saveRawLKr:saveRawCsI1:saveRawCsI2:saveRawMu:save_cal_rec:save_ktrrec");
    fw->modList();                                          //���������� ����� ���� ������� � �� ��������� (�������/��������)
    fw->beginJob();
    fw->proceed(5);                                         //����� �������������� �������, 0 - ���
    fw->endJob();
/*
    int nev,readerr, closeerr, openerr;
    int Nevents=500000;
    char fname[80];

    sprintf(fname,"/home/ovtin/simulation/simout/sim000001.dat");
    kedr_open_nat(fname,&openerr);
    if(openerr){
      printf("open error for %s\n",fname);
      exit(1);
    }
    printf("File %s is opened\n",fname);

    readerr=0;
    nev=0;

    while(nev<Nevents){
    kedr_read_nat(&readerr);
    if(readerr) break;

     nev++;
    }
    kedr_close_nat(&closeerr);
    printf("\nFile %s is closed\n",fname);
*/
    return 0;
}
