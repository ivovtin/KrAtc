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

#include "KrAtc/AtcRec.hh"
#include "KdATCCalib/AtcLed.hh"
#include "KrAtc/AtcHit.hh"
#include "KrAtc/AtcPar.hh"
#include "KrAtc/DCTrackExt.h"                    //описываются пересечения с цилиндром и плоскостями трека (прямой и спирали)
#include "KrAtc/atc_geometry.h"            //геометрия счетчиков

#include "KaFramework/mknatfilelist.h"

/*
#include "VDDCRec/kdcpar.h"
#include "VDDCRec/kvdpar.h"
#include "VDDCRec/ktracks.h"
#include "VDDCRec/VDDCRec.hh"
#include "VDDCRec/mtofhits.h"
#include "VDDCRec/kdcswitches.h"
#include "VDDCRec/kdcvd.h"
#include "VDDCRec/khits.h"
#include "VDDCRec/ToFTrack.hh"
extern "C" {
#include "VDDCRec/ktracks.h"
#include "VDDCRec/mtofhits.h"
}

#include "KEmcRec/emc_system.h"
#include "KEmcRec/emc_struct.h"

#include "KrVDDCMu/dcmu.h"
#include "KrMu/mu_system.h"
#include "KrMu/mu_event.h"

#include "KrToF/tof_system.h"

#include "KDisplay/kdisplay_event.h"
*/
//#include "TROOT.h"
//#include "TRint.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TGraphErrors.h"
//#include "TCanvas.h"


//  цепочка создается в конструкторе класса AppFramework с помощью его член-функции addModule
AppFramework::AppFramework(const char* sequence){
	addModule( new ReadNat     ("ReadNat", "Data reading"              ) );       //добавление модуля ReadNat
//        addModule( new VDDCRec     ("VDDCRec", "Tracks reconstruction"     ) );
	addModule( new AtcRec      ("AtcRec","ATC amplitude reconstruction") );
	reorder(sequence);             //функция переупорядочевания модулей
}

static const char* progname;
static const char* optstring="m:n:v";

void Usage(int status)
{
	cout<<"Usage: "<<progname<<" -m mode -n Nevents -v run|nat_file|run_list_file...\n"
		<<"Test ATC reconstruction\n"
		<<"Options:\n"
		<<"  -m mode     One of the AtcRec modes: test,run,led,gen (default: run)\n"
		<<"  -n Nevents  Number of events to process (default: all events)\n"
		<<"  -v          Be verbose\n"
		<<endl;
	exit(status);
}

int main(int argc, char* argv[])
{
     progname=argv[0];

	if( argc<2 ) Usage(0);        //использовать help

	int mode=AtcRec::RUN, dbglevel=1;                                   //по умолчанию значения
	long Nev=0;

//----------------- Process options -----------------//
	int opt;
	while( (opt=getopt(argc,argv,optstring))>0 ) {
		if( opt=='m' ) {
			string strmode=optarg;
			if     ( strmode=="test" ) mode=AtcRec::TEST;        //задаем режим - ген, лед, заход, либо тест для проверки
			else if( strmode=="run" )  mode=AtcRec::RUN;
			else if( strmode=="led" )  mode=AtcRec::LEDCAL;
			else if( strmode=="gen" )  mode=AtcRec::GENCAL;
			else {
				cerr<<"Unknown AtcRec mode "<<optarg<<endl;
				Usage(1);
			}
		} else if( opt=='n' ) {
			Nev=atoi(optarg);                                     //задаем число событий для обработки
			if( Nev<0 ) {
				cerr<<"Illegal number of events "<<Nev<<endl;
				return 1;
			}
		} else if( opt=='v' ) {                                      //включение verbose
		    dbglevel=5;                    //debugging parameter:
	                                           // 0 - do not print any messages
	                                           // 1 - print run and job statistics
	                                           // 2 - print event processing errors
	                                           // 3 - show a every event summary
	                                           // 4 - print hit info
                                                   // 5 - print everything, including fitting
		} else {
			Usage(1);
		}
	}

	flist_exev_t filelist;

	int nf=mknatfilelist(argc-optind,&argv[optind],filelist,1);    //считываем список nat файлов -  номер nat файла

	cout<<nf<<" files to process"<<endl;                        //число файлов для обработки
	if( !nf ) return -1;

	AppFramework* fw=new AppFramework("ReadNat->AtcRec");    //AppFramework- осущ-ет управление всеми модулями

	fw->verbose("AppFramework","cerr cout on");              //измен-е кол-ва выдаваемых модулем сообщений
	fw->verbose("AppFramework","clog off");
	fw->verbose("ReadNat","cerr cout clog on");
	fw->verbose("AtcRec","cerr cout clog on");

	fw->modList();                                          //напечатать имена всех модулей и их состояние (активен/пассивен)

	ostringstream sbuf;
	string empty;

	sbuf<<mode;                                            //режим - ген., свет, заход
	fw->modify("AtcRec::Mode",sbuf.str().c_str());         //изменить значение параметра модуля
	sbuf.str(empty);

	sbuf<<dbglevel;
	fw->modify("AtcRec::Debug",sbuf.str().c_str());        //изменить значение параметра модуля - verbose

	fw->modify("AtcRec::DoNotFit","0");
	fw->modify("AtcRec::MaxAmpDevSig","3.0");
	fw->modify("AtcRec::LoadCalibFromDB","1");
	fw->modify("AtcRec::MinDeltaT","53");
	fw->modify("AtcRec::MaxDeltaT","234");
	fw->modify("AtcRec::KdbReadMode","3");
	fw->modify("AtcRec::SkipDamagedEvent","0");


    if( Nev ) cout<<"Process "<<Nev<<" events"<<endl;         //вывод числа событий в заходе  - которое сами задаем для обр-ки

    long NevLeft=Nev;


	fw->beginJob();                                      //вызывается в начале работы для инициализации каждого модуля

	flist_exev_t::const_iterator iter=filelist.begin();     //первый нат файл

    ofstream list("out_file.txt",ios::out);                    //имя файла в который делаем вывод
    list<<"#hit.print()\n";

	while( iter!=filelist.end() ) {                       //обработка пока не закончатся nat файлы (список)
		fw->modify("ReadNat::input",(iter->first).c_str());
		long n=fw->proceed(NevLeft);                      //обработать NevLeft событий (если NevLeft=0 то обраб-ем все события). Возвращает число обработанных событий
		if( Nev ) {                                   //если есть события
			NevLeft-=n;
			if( NevLeft<=0 ) {
			    AppEvent temp;
			    fw->endRun(temp);                     //завершение захода
				break;
			}
		}

//	AtcHits *phits;
 /*       AppEvent event;
    if( event.get("AtcHits",phits) )
	return AppResult(AppResult::ERROR|AppResult::SKIP,"Can't get AtcHits");
*/
//    AtcHits::const_iterator first=phits->begin(), last=phits->end();

//========================================================================================================
      for(int i=1; i<160; i++)              //цикл по счетчикам
   {
     AtcHit hit(i);
    // AtcHit hit(35,a);        //вывод амп-ды конкретного сч-ка
    hit.rec(0);
    // hit.getType();                                 //тип события
   hit.print();                                  //вывод номера счетчика и амплитуд
   list<<"phits"<<endl;                        //то что выводим в файл
    }
//==========================================================================================================


		iter++;                                       //переходим на обработку следующего нат файла
     }

   	fw->endJob();                                         //вызывается в конце работы

    delete fw;

	return 0;
}
