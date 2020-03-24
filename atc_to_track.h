#ifndef ATC_TO_TRACK_H
# define ATC_TO_TRACK_H

# include "KcSys/fortran.h"
# include "AtcGlobals.h"

struct ATC_TRACK_INFO
{
	int track;                                                /*number of tracks (same as eTracksAll)*/
	int ncnt;                                                 /*number of counters crossed by the track*/
	int cnt[ATC_NCROS];                                       /*list of counter numbers on the track 1..160*/
	float npe[ATC_NCNT],npen[ATC_NCNT];         	          /*the hit amplitude in number of photoelectrons for each counter*/
	int wlshit[ATC_NCNT];                                     //WLS crossed
	int nearwls[ATC_NCNT];				          //track goes close to WLS
	float tlen[ATC_NCNT], pathwls[ATC_NCNT];		  //track length in aerogel volume and that in WLS

        //пересечение треком различных областей счетчика (0 или 1)
	int aerogel_REGION[ATC_NCNT];
	int aerogel_REGION0[ATC_NCNT];                  //ATC_DOUBLE_CROSS
	int aerogel_REGION5[ATC_NCNT];
	int aerogel_REGION20[ATC_NCNT];
	int active_REGION[ATC_NCNT];
	int active_REGION0[ATC_NCNT];
	int active_REGION5[ATC_NCNT];
	int active_REGION20[ATC_NCNT];
	int testreg[ATC_NCNT];				              //cosmic test region of counter crossed at total thickness
	int single_aerogel_REGION[ATC_NCNT];
	int single_aerogel_REGION0[ATC_NCNT];           //ATC_SINGLE_CROSS
	int single_aerogel_REGION5[ATC_NCNT];
	int single_aerogel_REGION20[ATC_NCNT];
	int single_active_REGION[ATC_NCNT];
	int single_active_REGION0[ATC_NCNT];
	int single_active_REGION5[ATC_NCNT];
	int single_active_REGION20[ATC_NCNT];
	int single_testreg[ATC_NCNT];				//cosmic test region of counter crossed at total thickness
};

extern struct ATC_TRACK_INFO atctrackinfo;

KC_SYS_FORTRAN_NAME(atctrackinfo, "atctrackinfo");          //common блоки в atcrec.inc

extern int atc_to_track(int t);

# ifdef __cplusplus
extern "C" int atc_info_(int* t);
# else
extern int atc_info_(int* t);
# endif


#endif
