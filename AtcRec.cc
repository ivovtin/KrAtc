#include <unistd.h>
#include <string>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "AppFramework/AppResult.hh"
#include "AppFramework/AppEvent.hh"
#include "ReadNat/re_def.h"
#include "ReadNat/rr_def.h"
#include "KDB/kdb.h"

#include "AtcRec.hh"

#include "AtcVFitter.h"
#include "A6PulseShape.h"
#include "AtcLeamaxFitter.h"
#ifndef ATC_NOROOT
# include "AtcMinuitFitter.h"
#endif

using namespace std;

const int AtcRec::chnl2addr[6][5] = {
{ 3, 0, 4, 8,12},
{ 7, 1, 5, 9,13},
{11, 2, 6,10,14},
{18,15,19,23,27},
{22,16,20,24,28},
{26,17,21,25,29} };

const char * const AtcRec::evTypeDesc[nEvTypes] = {
	" Empty","NoAddr","NoData","DupAdd","NotMap","  NoDT","IlAddr","IlData"," DTout","RecErr"
};

static int KedrEvent;

AtcRec::AtcRec(const char *nm, const char *descr) :
	AppModule(nm,descr),
	recConfFile(*this,"RecConfFile","/spool/users/atc/CalibPar/cntpar.atc"),
	dTcalFile(*this,"DTcalFile","dtrange.atc"),
	offCntFile(*this,"OffCntFile","offcnt.atc"),
	mode(*this,"Mode",RUN),
	digitMode(*this,"DigitMode",AtcHit::READ5),
	kdbReadMode(*this,"KdbReadMode",KDB_NEAREST),
	debug(*this,"Debug",0),
	skipDamagedEvent(*this,"SkipDamagedEvent",0),
	doNotFit(*this,"DoNotFit",0),
	useDeltaT(*this,"UseDeltaT",1),
	useTcutForCorrection(*this,"UseTcutForCorrection",0),
    loadCalibFromDB(*this,"LoadCalibFromDB",1),
    useLeamaxFitter(*this,"UseLeamaxFitter",0),
	printSkipCount(*this,"PrintSkipCount",100),
	badEventMask(*this,"BadEventMask",0),
	maxAmpDevSig(*this,"MaxAmpDevSig",3.0),
	minDeltaT(*this,"MinDeltaT",53.),
	maxDeltaT(*this,"MaxDeltaT",234.),
	spikeMinDev(*this,"SpikeMinDev",3.),
    spikeMinRatio(*this,"SpikeMinRatio",2.5),
    leamaxFitter(0), minuitFitter(0), fitter(0)
{}
AtcRec::~AtcRec()
{
	if( fitter ) delete fitter;
}

AppResult AtcRec::beginJob(AppEvent& event)
{
	cout<<"Job initialization"<<endl;

	int verboseFitting=(int)debug>=5?1:0;

#ifdef ATC_NOFIT
	doNotFit=true;
#endif

	leamaxFitter=new AtcLeamaxFitter(verboseFitting);    //LEAMAX - Constrained Non-Linear Least Squares and Maximum Likelihood Estimation
#ifndef ATC_NOROOT
	minuitFitter=new AtcMinuitFitter(verboseFitting);    //MINUIT - Function Minimization and Error Analysis
#endif //ATC_NOROOT

	if( (bool)doNotFit )
		cout<<"No pulse fitting"<<endl;
	else {
		if( (bool)useLeamaxFitter ) {
			cout<<"Leamax package is used for pulse fitting"<<endl;
			fitter=leamaxFitter;
		} else {
			cout<<"Minuit is used for pulse fitting"<<endl;
			fitter=minuitFitter;
		}

		leamaxFitter->releaseAll(); minuitFitter->releaseAll();
		leamaxFitter->fix(ATCPAR_P); minuitFitter->fix(ATCPAR_P);  //always fix P parameter at 5.0

		switch((int)mode) {
			case TEST:
				//test mode - fix nothing but P
				break;
			case RUN:
				//effect - fix pedestal only
				leamaxFitter->fix(ATCPAR_A0);
				if( minuitFitter ) minuitFitter->fix(ATCPAR_A0);
				break;
			case LEDCAL: case BGCAL:
				//LED or Background calibration - fix Pedestal & Tau
				leamaxFitter->fix(ATCPAR_A0);
				if( minuitFitter ) minuitFitter->fix(ATCPAR_A0);
				leamaxFitter->fix(ATCPAR_TAU);
				if( minuitFitter ) minuitFitter->fix(ATCPAR_TAU);
				break;
			case GENCAL:
				//Charge calibration - fix pedestal only
				leamaxFitter->fix(ATCPAR_A0);
				if( minuitFitter ) minuitFitter->fix(ATCPAR_A0);
				break;
			default:
				mode=RUN;
				leamaxFitter->fix(ATCPAR_A0);
				if( minuitFitter ) minuitFitter->fix(ATCPAR_A0);
				return AppResult(AppResult::WARNING,"Reconstruction mode is undefined. Choose default for effect.");
		};
	}

	lcnMap.clear();

	AtcHit::badEvent|=(int)badEventMask;
	AtcHit::setDigitMode((int)digitMode); //set default digit mode
	AtcHit::leamaxFitter=leamaxFitter;
	AtcHit::minuitFitter=minuitFitter;
	AtcHit::fitter=fitter; //set chosen fitter for AtcHit
	AtcHit::maxAmpDevSig=(float)maxAmpDevSig;
	AtcHit::spikeMinDev=(float)spikeMinDev;
	AtcHit::spikeMinRatio=(float)spikeMinRatio;
	AtcHit::verbose=(int)debug>=4?1:0;

	atcPar.minDeltaT=(float)minDeltaT;
	atcPar.maxDeltaT=(float)maxDeltaT;

	lastPedTime=0;
	last1peTime=0;
	lastGenTime=0;

	LoadRecConf();

	jobStat.clear();

	return AppResult();
}

AppResult AtcRec::endJob(AppEvent& event)
{
	if( (int)debug ) {
		cout<<"Job statistics:"<<endl;
		jobStat.print();
	}

	cout<<"End of job"<<endl;

	return AppResult();
}

AppResult AtcRec::beginRun(AppEvent& event)
{
	int m=(int)mode;
	bool norunnum=true;
	run=0;

	if( (norunnum=event.get("runNumber",run)) && m==RUN )
		return AppResult(AppResult::ERROR|AppResult::SKIP,"Can't get run number");

	cout<<"Run "<<run<<" initialization"<<endl;

	const RunRecord *runrec=0;
	if( m==LEDCAL || m==BGCAL || m==GENCAL )
		if( event.get("RunRecord",runrec) ) cerr<<"No run record"<<endl;

	//Try to determine digitizing mode automatically where possible
	if( m==RUN ) {
		if( run>4060 ) //4060 is the last run of the spring-summer season in 2004
			digitMode=AtcHit::READ5;
		else            //in fall of 2004 we went to the new A6 digitizing algorithm
			digitMode=AtcHit::READ4p1;
		AtcHit::setDigitMode((int)digitMode);
	} else if ( m==LEDCAL || m==BGCAL || m==GENCAL ) {
		if( run>=410 ) //410 - first number of calibration in fall 2004, though
			digitMode=AtcHit::READ5;   //there were no run record in nat file until December 2004
		else
			digitMode=AtcHit::READ4p1;
		AtcHit::setDigitMode((int)digitMode);

		if( runrec ) {
			struct tm btm;
			btm.tm_year=runrec->Header.Date.Year+100;
			btm.tm_mon=runrec->Header.Date.Month-1;
			btm.tm_mday=runrec->Header.Date.Day;
			btm.tm_hour=runrec->Header.Time.Hour;
			btm.tm_min=runrec->Header.Time.Minute;
			btm.tm_sec=0;
			runTime=mktime(&btm);
			cout<<"Begin time of calibration "<<run<<" of type "<<(int)runrec->Header.RunType<<": "
				<<ctime(&runTime)<<endl;
		}
	}

	if( (bool)loadCalibFromDB ) {
		if( LoadCalibFromDB() )
			return AppResult(AppResult::ERROR|AppResult::SKIP,"Failed to load calibration from DB");
	}

	runStat.clear();

	return AppResult();
}

AppResult AtcRec::endRun(AppEvent& event)
{
	if( (int)debug ) {
		cout<<"Run "<<run<<" statistics:"<<endl;
		runStat.print();
	}

	jobStat+=runStat;

	cout<<"End of run"<<endl;

	return AppResult();
}

AppResult AtcRec::event(AppEvent& event)
{
	static int rc;
	static pair<const short*,int> dir;

	runStat.total++;

	hits.clear();

	if( event.get("AtcRawRecord",dir) )
		return AppResult(AppResult::ERROR|AppResult::SKIP,"Can't get ATC record");
	if( event.get("eventNumber",KedrEvent) )
		return AppResult(AppResult::ERROR|AppResult::SKIP,"Can't get event number");

	//TODO: header parsing

	bool printErrors=(int)debug>=2?true:false;
	bool printBrief=(int)debug>=3?true:false;
	bool printAll=(int)debug>=4?true:false;

	if( printBrief ) {
		cout<<"Event "<<KedrEvent<<"\n";
		cout<<"ATC record of "<<dir.second<<" words"<<endl;
	}

	if( !dir.second ) {
		runStat.classified[EC_EMPTYREC]++;
		return AppResult();
	}

	addrEventRead.reset();

	const short *record = dir.first;
	const short *record_end = record + dir.second;

	int iamp, cnt, lcn=0, naddr=0, addr;
	short data;
	bool noaddress=false, nodata=false, nodeltat=false;
	bool illaddress=false, illdata=false;
	bool dupaddr=false, notmapped=false, bad=false;

//----------------- Read raw ATC record -----------------
	while( record < record_end ) {
		lcn = -(*record++);                        //номер логического канала
		data =  *record++;                         //данные

        addr = lcn-22001;                       //lcn - первый для АЧС 22001, последний 22841  (30(6*5)*28=840)

		if( printAll ) {
			::cout<<setw(5)<<lcn<<setw(5)<<data<<"  ";
			naddr++;                                        //число адресов
			if( naddr%6==0 ) ::cout<<endl;
		}
		if( lcn<0 ) {
			if( printErrors )
				cerr<<"ERROR Event "<<KedrEvent<<": Word "<<-lcn<<" is not an address"<<endl;
			noaddress=true;
			record--;
			continue;
		}
		if( addr<0 || addr>=2000 ) {
			if( printErrors )
				cerr<<"ERROR Event "<<KedrEvent<<": LCN="<<lcn<<" does not belong to ATC"<<endl;
			illaddress=true;
			break; //severe damage
		}
		if( data<0 ) {
			if( printErrors )
				cerr<<"ERROR Event "<<KedrEvent<<": Word "<<data<<" is not a data word"<<endl;
			nodata=true;
			record--;
			continue;
		}
		if( addrEventRead[addr] ) {
			if( printErrors )
				cerr<<"ERROR Event "<<KedrEvent<<": LCN "<<lcn<<" is already read for this event"<<endl;
			dupaddr=true;
			continue;
		}
		addrEventRead.set(addr);                               //установка прочитанного канала
		if( lcn==23000 ) //delta T data
		{
			if( data>255 ) {
				if( printErrors )
					cerr<<"ERROR Event "<<KedrEvent<<": Illegal value of DeltaT="<<data<<endl;
				illdata=true;
				continue;
			}
			hits.deltaT = data;
			continue;
		} else if( data>1023 ) {
			if( printErrors )
				cerr<<"ERROR Event "<<KedrEvent<<": Illegal value of amplitude="<<data<<endl;
			illdata=true;
			continue;
		}
		if( lcnMap.find(lcn)==lcnMap.end() ) { //may occur for earliest runs
			if( printErrors )
				cerr<<"ERROR Event "<<KedrEvent<<": LCN "<<lcn<<" not mapped to a counter"<<endl;
            notmapped=true;
			continue;
		}

		cnt=lcnMap[lcn]/5+1; // physics counter number
		iamp=lcnMap[lcn]%5;  // amplitude number 0-4

		if( hits.find(cnt)==hits.end() ) { //if this counter not read yet
			pair<int,AtcHit> curhit(cnt,AtcHit(cnt));
			curhit.second.a[iamp]=data;
			hits.insert(curhit); //add the hit to the collection of hits
		}
		else //if a hit of the counter already exists in the event
			hits[cnt].a[iamp]=data;
	} //while
	if( printAll && naddr%6!=0 ) ::cout<<endl;

	if( hits.deltaT==0 ) {
		nodeltat=true;
		runStat.classified[EC_NODELTAT]++;
		if( (bool)useDeltaT && printErrors )
			cerr<<"ERROR Event "<<KedrEvent<<": No DeltaT's address in record"<<endl;
	}

	if( notmapped )
		runStat.classified[EC_NOTMAPPED]++;

	if( illaddress || noaddress || nodata || dupaddr || illdata || nodeltat && (bool)useDeltaT ) {
		bad=true;
		if( noaddress )		runStat.classified[EC_NOADDR]++;
		if( nodata ) 		runStat.classified[EC_NODATA]++;
		if( dupaddr ) 		runStat.classified[EC_DUPADDR]++;
		if( illaddress )	runStat.classified[EC_ILLADDR]++;
		if( illdata )		runStat.classified[EC_ILLDATA]++;
		if( (bool)skipDamagedEvent ) {
			runStat.bad++;
			runStat.skipped++;
			if( printErrors )
				cerr<<"SKIP Event "<<KedrEvent<<": ATC record damaged"<<endl;
			return AppResult(AppResult::SKIP);
		}
		hits.ev_type=1; //if proceeded warn following modules that record is damaged
	}

//----------------- Reconstruction -----------------
	if( !nodeltat && (bool)useDeltaT ) {
		if( hits.deltaT<atcPar.minDeltaT-5 || hits.deltaT>atcPar.maxDeltaT+5 ) {
			runStat.classified[EC_DTOUT]++;
			hits.ev_type=2;
			if( printErrors ) {
				cerr<<"ERROR Event "<<KedrEvent<<": DeltaT="<<hits.deltaT<<" is out of range ("
					<<(atcPar.minDeltaT-5)<<", "<<(atcPar.maxDeltaT+5)<<")"<<endl;
			}
		}
		//make DeltaT correction in any case
		hits.deltaTfinal=(hits.deltaT-minDeltaT)/(maxDeltaT-minDeltaT)-0.5;
		if( printBrief )
			cout<<"DeltaT: Raw="<<hits.deltaT
				<<" Normalized="<<setprecision(3)<<hits.deltaTfinal<<endl;
	}

	if( printBrief )
		cout<<addrEventRead.count()<<" addresses read, "<<hits.size()<<" hits found"<<endl;

	bool errinrec=false;
	AtcHits::iterator cur=hits.begin(), last=hits.end();
	int nbad=0, nest=0, nfit=0;

	for( ; cur!=last; cur++)        //цикл по хитам
	{
		cnt=cur->first;
   		AtcHit *ahit=&cur->second;

		rc=ahit->rec(doNotFit);

		if( rc<0 ) //initializing error
			return AppResult(AppResult::STOP|AppResult::ERROR,"Problem with initialization");

		if( rc==1 ) { //hit corrupted
			nbad++;
			if( printErrors )
				cerr<<"ERROR Event "<<KedrEvent<<": Bad amplitude in counter "<<cnt<<endl;
			errinrec=true;
			//do not proceed further if the event will be skipped
			if( (bool)skipDamagedEvent ) break;
			continue;
		}

		//make timing correction if the hit is at least estimated
		if( ahit->checkType(AtcHit::Estimated) ) {
			nest++;
			if( ahit->checkType(AtcHit::Fitted) ) nfit++;

			ahit->Tfinal=ahit->T+hits.deltaTfinal;

			if( useTcutForCorrection )
			{
				if( cntPar[cnt].tCutRight!=0. &&
					cntPar[cnt].dtCutRight!=0. &&
					ahit->Tfinal > cntPar[cnt].tCutRight &&
					hits.deltaTfinal > cntPar[cnt].dtCutRight ) ahit->Tfinal-=1.;
				if( cntPar[cnt].tCutLeft!=0. &&
					cntPar[cnt].dtCutLeft!=0. &&
					ahit->Tfinal < cntPar[cnt].tCutLeft &&
					hits.deltaTfinal < cntPar[cnt].dtCutLeft ) ahit->Tfinal+=1.;
			}
		}
		if( printAll ) ahit->print();
	} //while

	if( printBrief )
		cout<<hits.size()<<" hits: "<<nbad<<" bad, "<<nest<<" have peak, "<<nfit<<" fitted"<<endl;

	if( errinrec ) {
		bad=true;
		runStat.classified[EC_ERRINREC]++;
		if( (bool)skipDamagedEvent ) {
			runStat.skipped++;
			if( printErrors )
				cerr<<"SKIP Event "<<KedrEvent<<": Hit(s) corrupted"<<endl;
			return AppResult(AppResult::SKIP);
		}
		hits.ev_type=1;
	}

	if( bad )
		runStat.bad++;

	if( event.put("AtcHits", &hits) )
		return AppResult(AppResult::ERROR|AppResult::SKIP,"Can't put reconstructed ATC data");

	runStat.processed++;

	return AppResult();
}
//загрузка крнфигурации реконструкции
int AtcRec::LoadRecConf()
{
	int i, c, ic, a6in, lcn;        //i - точка срабатывания(5), ic - номре счетчика, lcn-Logical Channel Number
//----------------- Reading reconstruction configuration -----------------
	cout<<"Loading ATC reconstruction parameters from "<<(string)recConfFile<<endl;
	ifstream in(recConfFile.c_str());
	if( !in.is_open() ) {
		cerr<<"Can not open "<<(string)recConfFile<<". Using default LCN mapping."<<endl;
		for(ic=0; ic<NATC; ic++) {
			for(i=0; i<5; i++) {
				a6in = ic+(ic/80)*4;
				lcn=22001+(a6in/6)*30+chnl2addr[a6in%6][i];          //a6in/6 - номер позиции платы А6 в крейте
				lcnMap[lcn]=ic*5+i; //store lcn to amplitude index mapping
			}
		}
	} else {
		int cnt,ncard,ch;
		float ped,dped,tau,dtau,p,dp;
		while(1) {
			//read a string:
			//c     - physics counter number 1-160, negative if counter is off
			//ncard - logical card number 1-28 as determined in ATC.KDS
			//ch    - channel number in the card
			//ped,dped - pedestal value and rms previously calibrated
			//tau,dtau - time constant parameter value and rms previously calibrated
			//p,dp  - power parameter value and rms previously calibrated
			in>>c>>ncard>>ch>>ped>>dped>>tau>>dtau>>p>>dp;

			if( in.eof() ) break;
			if( in.fail() ) {
				cerr<<"ERROR Reading "<<recConfFile.c_str()<<" failed\n";
				return 2;
			}

			cnt=abs(c);
			if( cnt>NATC || cnt<1 ) {
				cerr<<"ERROR Wrong counter number cnt="<<cnt<<endl;
				continue;
			}
			if( ncard<1 || ch<1 || ch>6 ) {
				cerr<<"ERROR Wrong A6 card geometry: cnt="<<cnt
					<<" NA6="<<ncard<<" Channel="<<ch<<endl;
				continue;
			}

			cntPar[cnt].cnt=cnt;
			if( c<0 ) cntPar[cnt].off=true;

			if( ped>=0. && ped<100. ) {
				if( dped>1 )
					cntPar[cnt].setfA0(ped,dped);
				else
					cntPar[cnt].setfA0(ped);
			}

			if( tau>0. )
				cntPar[cnt].setfTau(tau,dtau);
			else
				cntPar[cnt].setfTau(1.0,0.1);

			//always fix P parameter at 5.0
			cntPar[cnt].setfP(5.0,0.0);

			for(i=0; i<5; i++) {
				lcn=22001+(ncard-1)*30+chnl2addr[ch-1][i];
				if(lcn<24000)
					lcnMap[lcn]=(cnt-1)*5+i;
				else
					cerr<<"ERROR Illegal LCN="<<lcn<<": cnt="<<cnt<<" NA6="<<ncard<<" Channel="<<ch<<endl;
			}
		}//while
	}
	in.close();
	in.clear(); // now one can reuse this stream again with another file

//----------------- Reading DeltaT range -----------------
	int r1, r2;
	in.open(dTcalFile.c_str());
	if( in.is_open() ) {
		cout<<"Reading DeltaT range from "<<(string)dTcalFile<<endl;
		in>>r1>>r2;
		if( in.fail() ) {
			cerr<<"ERROR Reading "<<dTcalFile<<" failed"<<endl;
		} else {
			if( r1==0 && r2==0 ) {
				useDeltaT=false;
				cout<<"Set DeltaT correction off"<<endl;
			} else if( r1<0 || r1>255 || r2<0 || r2>255 || r1>r2 ) {
				cerr<<"ERROR Wrong DeltaT range: ("<<r1<<","<<r2<<")"<<endl;
			} else {
				atcPar.minDeltaT=r1;
				atcPar.maxDeltaT=r2;
				cout<<"DeltaT range read: ("<<r1<<","<<r2<<")"<<endl;
			}
		}
	}
	in.close();
	in.clear();

//----------------- Reading off-counters list -----------------
	in.open(offCntFile.c_str());
	if( in.is_open() ) {
		cout<<"Loading off-counters list from "<<(string)offCntFile<<endl;
		while(1) {
			in>>c;
			if( in.eof() ) break;
			if( in.fail() ) {
				cerr<<"ERROR Reading "<<offCntFile<<" failed"<<endl;
				break;
			}
			if( c<1 || c>NATC ) {
				cerr<<"Wrong counter number cnt="<<c<<endl;
				continue;
			}
			cntPar[c].off=true;
		}
	}
	in.close();
	in.clear();

//----------------- Reading time cut ranges -----------------
	if( !(bool)useTcutForCorrection ) return 0;
	cout<<"Loading time cut parameters from "<<(string)tCutParFile<<endl;
	in.open(tCutParFile.c_str());
	if( !in.is_open() ) {
		cerr<<"ERROR Can not open "<<(string)tCutParFile<<endl;
		useTcutForCorrection=false;
		return 1;
	}
	float tl,dtl,tr,dtr;
	while(1) {
		//read a string:
		//c     - physics counter number 1-160
		//tl    - time cut left
		//dtl   - deltaT cut left
		//tr    - time cut right
		//dtr   - deltaT cut right
		in>>c>>tl>>dtl>>tr>>dtr;
		if( in.eof() ) break;
		if( in.fail() ) {
			cerr<<"ERROR Reading "<<tCutParFile<<" failed"<<endl;
			useTcutForCorrection=false;
			return 2;
		}
		if( c<1 || c>NATC ) {
			cerr<<"ERROR Wrong counter number cnt="<<c<<endl;
			continue;
		}
		cntPar[c].setTdTcut(tl,dtl,tr,dtr);
	}
	in.close();
	in.clear();

	return 0;
}
//загрузка калибровки из БД
int AtcRec::LoadCalibFromDB()
{
	static int pedCalN=0, ledCalN=0, genCalN=0;

	int pedDB[Ped_len], ledDB[LED_len], genDB[Gen_len], dtDB [DT_len];
	bool usePrevPed=false, usePrev1pe=false, usePrevTau=false;
	time_t caltime;

	noPedInDB=true;
	no1peInDB=true;
	noTauInDB=true;
	noDTInDB =true;

	KDBconn* connection=kdb_open();
	if( !connection ) {
		cerr<<"Can not establish connection to database"<<endl;
		return -1;
	}

	switch( (int)kdbReadMode ) {
		case 1: cout<<"KDB read mode: KDB_VALID"<<endl;
			kdb_set_readmode(connection,KDB_VALID);
			break;
		case 2: cout<<"KDB read mode: KDB_NEAREST"<<endl;
			kdb_set_readmode(connection,KDB_NEAREST);
			break;
		case 3: cout<<"KDB read mode: KDB_NEXT"<<endl;
			kdb_set_readmode(connection,KDB_NEXT);
			break;
		case 4: cout<<"KDB read mode: KDB_PREVIOUS"<<endl;
			kdb_set_readmode(connection,KDB_PREVIOUS);
			break;
		default: cout<<"KDB read mode "<<(int)kdbReadMode<<" is unknown. Use KDB_VALID."<<endl;
			kdb_set_readmode(connection,KDB_VALID);
	}

	int m=(int)mode;

	if( m==RUN ) {
		runTime=kdb_run_get_begin_time(connection, run);
		cout<<"Begin time of "<<run<<": "<<timestamp(runTime)<<endl;
	}

	if( runTime<=0 ) {
		kdb_close(connection);
		return m==TEST?0:1;
	}

	//load calibration only if it is different from previous one
	caltime=kdb_get_begin_time(connection,LED_ID,runTime);
	if( last1peTime==0 || last1peTime!=caltime ) {
		if( kdb_read(connection, LED_ID, runTime, ledDB, LED_len) ) {
			last1peTime=caltime;
			ledCalN=ledDB[LED_len-1];
			no1peInDB=false;
			cout<<"1pe calibration "<<ledCalN<<" loaded ("<<timestamp(caltime)<<")"<<endl;
		} else {
			cerr<<"Failed to get 1pe calibration"<<endl;
		}
	} else {
		usePrev1pe=true;
		no1peInDB=false;
		cout<<"Use loaded 1pe calibration "<<ledCalN<<" ("<<timestamp(caltime)<<")"<<endl;
	}

	caltime=kdb_get_begin_time(connection,Ped_ID,runTime);
	if( lastPedTime==0 || lastPedTime!=caltime ) {
		if( kdb_read(connection, Ped_ID, runTime, pedDB, Ped_len) ) {
			pedCalN=pedDB[1];
			lastPedTime=caltime;
			noPedInDB=false;
			cout<<"Pedestal calibration "<<pedCalN<<" loaded ("<<timestamp(caltime)<<")"<<endl;
		} else {
			cerr<<"Failed to get pedestals calibration"<<endl;
		}
	} else {
		usePrevPed=true;
		noPedInDB=false;
		cout<<"Use loaded pedestal calibration "<<pedCalN<<" ("<<timestamp(caltime)<<")"<<endl;
	}

	if( m!=GENCAL ) {
		caltime=kdb_get_begin_time(connection,Gen_ID,runTime);
		if( lastGenTime==0 || lastGenTime!=caltime ) {
			if( kdb_read(connection, Gen_ID, runTime, genDB, Gen_len) ) {
				genCalN=genDB[Gen_len-1];
				lastGenTime=caltime;
				noTauInDB=false;
				cout<<"Charge calibration "<<genCalN<<" loaded ("<<timestamp(caltime)<<")"<<endl;
			} else {
				cerr<<"Failed to get tau calibration"<<endl;
			}
		} else {
			usePrevTau=true;
			noTauInDB=false;
			cout<<"Use loaded charge calibration "<<genCalN<<" ("<<timestamp(caltime)<<")"<<endl;
		}
	}

	if( kdb_read(connection, DT_ID, runTime, dtDB, DT_len) ) {
		noDTInDB=false;
		cout<<"DeltaT calibration: "<<dtDB[0]<<" "<<dtDB[1]<<" ("<<timestamp(caltime)<<")"<<endl;
	} else {
		cerr<<"Failed to get DeltaT calibration"<<endl;
	}

	kdb_close(connection);

	if( !noDTInDB ) {
		if( dtDB[0]>0 && dtDB[0]<255 && dtDB[1]>0 && dtDB[1]<=255 && dtDB[0]<dtDB[1] ) {
			atcPar.minDeltaT=dtDB[0];
			atcPar.maxDeltaT=dtDB[1];
		} else {
			cerr<<"DeltaT calibration loaded is invalid and was ignored"<<endl;
		}
	}

	if( !noPedInDB && !usePrevPed ) {
		for(int i=0, c=1; i<NATC; i++, c++) {
			float ped=pedDB[i*2+2]/1000.;
			float dped=pedDB[i*2+3]/1000.;
			if( ped>0. && ped<100. ) {
				if( dped>1 )
					cntPar[c].setfA0(ped,dped);
				else
					cntPar[c].setfA0(ped);
			}
		}
	}

	if( !noTauInDB && !usePrevTau ) {
		for(int i=0, c=1; i<NATC; i++, c++) {
			float tau=genDB[i*8]/100.;
			float dtau=genDB[i*8+1]/100.;
			if( tau>0. ) cntPar[c].setfTau(tau,dtau>=0.?dtau:0.);
		}
	}

	if( !no1peInDB && !usePrev1pe ) {
		int nbad=0;
		for(int i=0, c=1; i<NATC; i++, c++) {
			float a1=ledDB[i*8]/100.;
			float da1=ledDB[i*8+1]/100.;
			if( a1>0. ) cntPar[c].setA1pe(a1,da1>=0.?da1:0.);
			else nbad++;
		}
		cout<<nbad<<" bad counters in database"<<endl;
	}

	//determine that loading failed
	if( m!=TEST && noPedInDB || m==RUN && no1peInDB ) return 1;

	return 0;
}

void AtcRec::EventStat::print(ostream& out) const
{
	out<<"Total number of events     = "<<setw(6)<<total<<"\n";
	out<<"Number of processed events = "<<setw(6)<<processed<<"  "<<setprecision(3)<<100.*processed/total<<"%"<<"\n";
	out<<"Number of bad events       = "<<setw(6)<<bad      <<"  "<<setprecision(3)<<100.*bad/total      <<"%"<<"\n";
	out<<"Number of skipped events   = "<<setw(6)<<skipped  <<"  "<<setprecision(3)<<100.*skipped/total  <<"%"<<"\n";
	out<<"Number of classified errors:\n";

	int i;
	for(i=0; i<nEvTypes; i++)
		out<<" "<<evTypeDesc[i]<<" ";
	out<<"\n";
	for(i=0; i<nEvTypes; i++)
		out<<setw(7)<<classified[i]<<" ";
	out<<"\n";
	for(i=0; i<nEvTypes; i++)
		out<<fixed<<setprecision(5)<<setw(7)<<100.*classified[i]/total<<" ";
	out<<endl;

}

