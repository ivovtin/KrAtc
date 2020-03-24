#ifndef ATC_REGIONS_H
# define ATC_REGIONS_H

# include "KcSys/fortran.h"
enum {
	ATC_AEROGEL_REGION,   //Aerogel area with offset from walls and shifter
	ATC_ACTIVE_REGION,    //Active area -> shifter+aerogel with offset from walls and shifter
	ATC_CRT_REGION,       //area trigger on cosmic tests
	ATC_AEROGEL_REGION0,  //different offset - 0, 5, 20 mm
	ATC_AEROGEL_REGION5,
	ATC_AEROGEL_REGION20,
	ATC_ACTIVE_REGION0,
	ATC_ACTIVE_REGION5,
	ATC_ACTIVE_REGION20
};
enum {
	ATC_DOUBLE_CROSS,     //different cases crosses counter
	ATC_SINGLE_CROSS,
	ATC_IN_CROSS,
	ATC_OUT_CROSS
};

extern int atc_track_cross_region(int t, int iatc,
	int regtype=ATC_AEROGEL_REGION, int crosstype=ATC_DOUBLE_CROSS);

# ifdef __cplusplus
extern "C" int atc_trreg_(int* nt, int* natc, int *regtype, int *crosstype);
# else
extern int atc_trreg_(int* nt, int* natc, int *regtype, int *crosstype);
# endif


#endif //ATC_REGIONS_H
