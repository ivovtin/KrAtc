#ifndef ATCREC_H
# define ATCREC_H

# include "KcSys/fortran.h"
# include "AtcGlobals.h"

/*Data for ATC counters bound to DC tracks*/
struct AtcTrackData {
	/*number of tracks (same as eTracksAll)*/
	int ntracks;
	/*number of counters crossed by the track*/
	unsigned char ncnt_on_track[ATC_NTRK];
	/*list of counter numbers on the track 1..160*/
        int cnt_cross[ATC_NTRK][ATC_NCROS];
        /*whether DC track is wrongly reconstructed as determined by ATC*/
	unsigned char illtrack[ATC_NTRK];
	/*whether WLS is crossed by the track in the counter*/
	unsigned char wls_hit[ATC_NTRK][ATC_NCROS];
	/*whether the track had possibility to cross WLS*/
	unsigned char near_wls[ATC_NTRK][ATC_NCROS];
	/*which end of the track crossed the counter: 0 or 1 (for cosmic run only)*/
	unsigned char track_end[ATC_NTRK][ATC_NCROS];
	/*track length in aerogel and WLS for the counter, mm*/
	float tlen_in_aer[ATC_NTRK][ATC_NCROS], tlen_in_wls[ATC_NTRK][ATC_NCROS];
	/*cylindrical rho coordinate of inner and outer crossings of the counter, mm*/
	float rin[ATC_NTRK][ATC_NCROS], rout[ATC_NTRK][ATC_NCROS];
	/*cylindrical phi coordinate of inner and outer crossings of the counter, rad*/
	float phiin[ATC_NTRK][ATC_NCROS], phiout[ATC_NTRK][ATC_NCROS];
	/*cylindrical z coordinate of inner and outer crossings of the counter, mm*/
	float zin[ATC_NTRK][ATC_NCROS], zout[ATC_NTRK][ATC_NCROS];
 	/*track phase for inner and outer crossings of the counter, rad*/
	float ph_in[ATC_NTRK][ATC_NCROS], ph_out[ATC_NTRK][ATC_NCROS];
	/*global coordinates of crossing WLS, mm&rad*/
	float rwls[ATC_NTRK][ATC_NCROS], phiwls[ATC_NTRK][ATC_NCROS], zwls[ATC_NTRK][ATC_NCROS];
	/*global coordinates of crossing barrel layers, mm&rad*/
	float rB[ATC_NTRK][8],phiB[ATC_NTRK][8],zB[ATC_NTRK][8];
	/*global coordinates of crossing endcap layers, mm&rad*/
	float rE[ATC_NTRK][8],phiE[ATC_NTRK][8],zE[ATC_NTRK][8];

	/*link counters to tracks, tracks indexed starting from 1*/
	unsigned char cnt_ntrk[ATC_NCNT], cnt_tracks[ATC_NCNT][ATC_NTRK];
};

#define TrackAtcData AtcTrackData

/*Data for each ATC counter.*/
struct AtcRecData {
	/*the hit amplitude in number of photoelectrons for each counter*/
	float npe[ATC_NCNT], dnpe[ATC_NCNT];
	/*time of the hit (DeltaT corrected)*/
	float time[ATC_NCNT];
	/*raw amplitude and time of the hit*/
	float amp[ATC_NCNT], rawtime[ATC_NCNT];
	/*Chi-square of pulse fit*/
	float chi2[ATC_NCNT];
	/*whether the counter is present in the event*/
	unsigned char trg[ATC_NCNT];
	/*hit type determined by AtcHit*/
	int htyp[ATC_NCNT];
	/*event corruption type, 0 if no damage seen*/
	int eventdamage;
	/*DeltaT reading 0..255*/
	int rawDeltaT;
	/*DeltaT correction in ticks -0.5..0.5*/
	float deltaT;
};

/*ATC calibration data. Mostly copied from AtcPar for exporting to Fortran.*/
struct AtcRunData {
	/*bad counter signatures*/
	int bad[ATC_NCNT];
	/*amplitude measurement errors*/
	float amp_err[ATC_NCNT];
	/*array of channel pedestals*/
	float a0[ATC_NCNT];
	/*array of tau parameters*/
	float tau[ATC_NCNT];
	/*array of p parameters*/
	float p[ATC_NCNT];
	/*one photoelectron amplitudes from LED calibration and errors*/
	float a1pe[ATC_NCNT], da1pe[ATC_NCNT];
	/*DeltaT calibration constants*/
	float minDeltaT, maxDeltaT;
	/*transversal and longitudinal coordinate resolution at ATC radius in mm.
	  used for determining proximity to WLS*/
	float dcSigmaT, dcSigmaZ;
};

/*tracks not to consider*/
extern unsigned char atc_droptr[ATC_NTRK];

/*if set to non-zero value do tracking to ATC*/
extern int atc_tracking;

# ifdef __cplusplus
class AtcRec;
extern AtcRec *atcRec;
# endif

extern struct AtcTrackData atc_track;
extern struct AtcRecData atc_rec;
extern struct AtcRunData atc_data;

KC_SYS_FORTRAN_NAME(atc_track,   "atc_track");
KC_SYS_FORTRAN_NAME(atc_rec,     "atc_rec");
KC_SYS_FORTRAN_NAME(atc_data,    "atc_data");
KC_SYS_FORTRAN_NAME(atc_droptr,  "atc_droptr");
KC_SYS_FORTRAN_NAME(atc_tracking,"atc_tracking");

# ifdef __cplusplus
extern "C" {
# endif
extern int atc_init();
extern int atc_run(int run);
extern int atc_event();
extern void atc_stop();
extern void rec_sim_event();
extern int get_atccalib_range_(int* range);
# ifdef __cplusplus
}
# endif

KC_SYS_FORTRAN_NAME(atc_init,  "atc_init");
KC_SYS_FORTRAN_NAME(atc_run,   "atc_run");
KC_SYS_FORTRAN_NAME(atc_event, "atc_event");
KC_SYS_FORTRAN_NAME(atc_stop,  "atc_stop");

#endif
