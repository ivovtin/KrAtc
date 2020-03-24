#ifndef AtcRec_hh
#define AtcRec_hh

#include <string>
#include <cmath>
#include <list>
#include <bitset>
#include <iosfwd>
#include <ctime>
#include "AppFramework/AppModule.hh"
#include "AppFramework/AppParameter.hh"
#include "KDB/kdb.h"

#include "AtcPar.hh"
#include "AtcHit.hh"

class AppResult;
class AppEvent;
class AtcVFitter;

class AtcRec : public AppModule        //AtcRec - класс, наследующий базовый класс  AppModule
{
	// Make AtcRec friendly to the C/Fortran interface
	friend int atc_init();                                //функция инициализации рек-ии
	friend int atc_run(int);                              //инициализация для захода
	friend int atc_event();                               //реконструкция одного события
	friend void atc_stop();                               //освобождение используемой памяти по окончании работы

public:
	// Parameters of reconstruction
	AppParameter<string> recConfFile;          //LCN mapping and rec parameters
	AppParameter<string> dTcalFile;            //DeltaT calibration file
	AppParameter<string> offCntFile;           //disabled counters
	AppParameter<string> tCutParFile;          //time cut ranges source
	AppParameter<int>    mode;                 //reconstruction mode:
	                                           // 0 - preliminary
	                                           // 1 - effect
	                                           // 2 - LED calibration
	                                           // 3 - GEN calibration
	                                           // 4 - Background calibration
	AppParameter<int>    digitMode;            //digitizing mode of A6: 0 - old (4+1), 1 - new (5)
	AppParameter<int>    kdbReadMode;          //reading mode of KDB
	                                           // 1 - nearest earlier (valid) calibration (KDB_VALID)
	                                           // 2 - nearest (earlier or later) calibration by time (KDB_NEAREST)
	                                           // 3 - next calibration excluding exactly matched (KDB_NEXT)
	                                           // 4 - previous calibration excluding exactly matched (KDB_PREVIOUS)
	AppParameter<int>    debug;                //debugging parameter:
	                                           // 0 - do not print any messages
	                                           // 1 - print run and job statistics
	                                           // 2 - print event processing errors
	                                           // 3 - show a every event summary
	                                           // 4 - print hit info
                                               // 5 - print everything, including fitting
	AppParameter<bool>   skipDamagedEvent;     //refuse events with damaged ATC record
	AppParameter<bool>   doNotFit;             //do not use fit, only estimate
	AppParameter<bool>   useDeltaT;            // use DeltaT for precise timing (default)
	AppParameter<bool>   useTcutForCorrection; // use (T,deltaT) cut for T correction
	AppParameter<bool>   loadCalibFromDB;      //whether to load calibration data from DB
	AppParameter<bool>   useLeamaxFitter;      //whether to use LEAMAX from CERNLIB package to fit pulses
	AppParameter<int>    printSkipCount;       //print every skip count event (obsolete, to be removed soon)
	AppParameter<int>    badEventMask;         //Mask of event types not fitted
	AppParameter<float>  maxAmpDevSig;         //Acceptable deviation of amplitude in sigma units
	AppParameter<float>  minDeltaT;            //Min deltaT value
	AppParameter<float>  maxDeltaT;            //Max deltaT value
	AppParameter<float>  spikeMinDev;          //Minimum peak deviation from substrate
	                                           // in number of errors
	AppParameter<float>  spikeMinRatio;        //Minimum spike amplitude to substrate span ratio

private:
	// Internal constants and variables

	static const int chnl2addr[6][5]; //A6 card channel address table

	hash_map<int,int> lcnMap;         //LCN map on counters

	AtcVFitter *leamaxFitter, *minuitFitter; //pointers to hit fitting managers
	AtcVFitter *fitter;                     //pointer to chosen hit fitting manager

	std::bitset<2000> addrEventRead; //bit set of addresses present in an event record

	int run;        //current run number
	time_t runTime; //start time of the current run
	time_t lastPedTime, //last begintime of pedestal calibration
	       last1peTime, //last begintime of 1pe calibration
	       lastGenTime; //last begintime of charge calibration
	bool noPedInDB, //no pedestal in database for current run
		 no1peInDB, //no 1pe calibration in database
	     noTauInDB, //no tau calibration in database
	     noDTInDB;  //no DeltaT calibration in database

public:
	// Hits container of current event
	AtcHits hits;

	// Supported reconstruction mode numbers
	enum RecModes { TEST, RUN, LEDCAL, GENCAL, BGCAL };

	// Event statistics

	static const int nEvTypes=10;
	static const char * const evTypeDesc[nEvTypes];
	enum { EC_EMPTYREC, EC_NOADDR, EC_NODATA, EC_DUPADDR, EC_NOTMAPPED,         //типы событий
	       EC_NODELTAT, EC_ILLADDR, EC_ILLDATA, EC_DTOUT, EC_ERRINREC };

	struct EventStat {                                                          //структура статистики
		long total, processed, bad, skipped;
		long classified[nEvTypes];

		EventStat& operator+=(const EventStat& r) {
			total+=r.total; processed+=r.processed; bad+=r.bad; skipped+=r.skipped;
			for(int i=0; i<nEvTypes; i++) classified[i]+=r.classified[i];
			return *this;
		}

		EventStat() { clear(); }

		void clear() {
			total=0; processed=0; bad=0; skipped=0;
			for(int i=0; i<nEvTypes; i++) classified[i]=0;
		}

		void print(std::ostream& out=std::cout) const;
	};

	EventStat jobStat, runStat;

    // Database tables info
	enum DBtableID     { Ped_ID=1302,  Gen_ID=1305,  LED_ID=1307,  DT_ID=1308 };
	enum DBtableLength { Ped_len=322,  Gen_len=1284, LED_len=1283, DT_len=2   };

private:
	// AppFramework called methods
	AppResult beginJob(AppEvent& event);     //для инициализации модуля
	AppResult beginRun(AppEvent& event);     //переход к след заходу
	AppResult event   (AppEvent& event);
	AppResult endRun  (AppEvent& event);
	AppResult endJob  (AppEvent& event);     //выз-ся в конце работы

private:
	// Class private functions
	int LoadRecConf();      //load reconstruction configuration file
	int LoadCalibFromDB();  //load calibration data from DB

	char* timestamp(time_t t); //return string representing given time value

public:
	AtcRec(const char *nm, const char *descr);
	~AtcRec();
};

inline char* AtcRec::timestamp(time_t t)
{
	static const char* time_fmt="%b %d %Y %H:%M:%S";
	static char strtime[50];

	struct tm* brtm=localtime(&t);
	strftime(strtime,50,time_fmt,brtm);

	return strtime;
}
#endif
