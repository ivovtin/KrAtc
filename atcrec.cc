#include "atcrec.h"

#include <iostream>
#include <string.h>
using std::cout;
using std::cerr;
using std::clog;
using std::endl;

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <numeric>

#include "AppFramework/AppResult.hh"
#include "AppFramework/AppEvent.hh"
#include "ReadNat/re_def.h"
#include "ReadNat/ss_def.h"
#include "VDDCRec/ktracks.h"
#include "VDDCRec/kdcpar.h"
#include "VDDCRec/kglobparam.h"
#include "AtcRec.hh"
#include "DCTrackExt.h"
#include "atc_geometry.h"
#include "atc_to_track.h"

#include "atc_sim.h"
#include "cernlib.h"
#include "ReadNat/mc_def.h"
#include "ReadNat/rr_def.h"
#include "ReadNat/mcatc_def.h"
#include "ReadNat/atcrawhitspar.h"
#include "ReadNat/read_nat_c.h"
#include "KcSys/fortran.h"
#include "/cern/pro/include/geant321/pilot.h"
#include "KDB/kdb.h"
#include <ctime>

//sim
extern mcatcEvent kedrmcatchits_;   //ReadNat
extern kscsim0_type kscsim0_;
extern void atc_sim_zero_raw_hits();

// Global variables section
AtcRec *atcRec=0;
struct AtcTrackData atc_track;
struct AtcRecData atc_rec;
struct AtcRunData atc_data;
struct ATC_TRACK_INFO atctrackinfo;
int atc_tracking=1;
unsigned char atc_droptr[ATC_NCROS];

// Local variables section
static int fMagField=-1;   //magnetic field flag (undefined first)
static AppResult result;
static float Tmin, Tmax;   //range of pulse times
static int runNumber=0;    //run number processed by atcrec

static int beam_rec(long eventNumber);
static int cosmic_rec(long eventNumber);

//pointer to event reconstruction function
static int (*event_rec)(long) = beam_rec;
double dx1, dy1, dz1;

//fill atc_data structure
static void fill_atc_data()
{
	int i, c;
	for(i=0, c=1; i<NATC; i++, c++) {
		atc_data.bad[i]=cntPar[c].isOffline()?1:0;
		atc_data.amp_err[i]=cntPar[c].getError();
		atc_data.a0[i]=cntPar[c].getfA0();
		atc_data.tau[i]=cntPar[c].getfTau();
		atc_data.p[i]=cntPar[c].getfP();
		atc_data.a1pe[i]=cntPar[c].getA1pe();
		atc_data.da1pe[i]=cntPar[c].getA1peErr();
	}
	atc_data.minDeltaT=atcPar.minDeltaT;
	atc_data.maxDeltaT=atcPar.maxDeltaT;
	//transversal and longitudinal DC resolution at ATC radius, mm
	atc_data.dcSigmaT=0.3;
	atc_data.dcSigmaZ=3;
}

int atc_init()
{
	if( atcRec ) return -1; //ATC reconstruction already inited

	atcRec = new AtcRec("AtcRec","ATC amplitude reconstruction");

	atcRec->mode=AtcRec::RUN;
	atcRec->badEventMask=0;
	atcRec->doNotFit=false;
	atcRec->skipDamagedEvent=true;
	atcRec->maxAmpDevSig=3.0;
	atcRec->loadCalibFromDB=true;

	AppEvent event;
	result=atcRec->beginJob(event);

	if( result.status&(AppResult::ERROR|AppResult::WARNING) )
		cerr<<result.message<<endl;
	if( result.status&AppResult::LOG )
		clog<<result.message<<endl;

	if( result.status&AppResult::STOP ) return -1;

	fill_atc_data();
	atc_fill_geometry();

	return 0;
}
void atc_stop()
{
	delete atcRec;
}
int atc_run(int run)
{
	if( run==runNumber ) return -1; //already inited for this run

	runNumber=run;

	AppEvent event;
	event.put("runNumber",run);

	result=atcRec->beginRun(event);

	if( (int)atcRec->digitMode==AtcHit::READ4p1 ) Tmin=0.0, Tmax=6.0;
	else Tmin=0.0, Tmax=4.0;

	if( result.status&AppResult::SKIP ) {
		cerr<<"ATC_RUN: Skipping run "<<run<<endl;
		return 1;
	}

      	fill_atc_data();

	return 0;
}

enum RC_CODES {
	RC_DCREC_ERROR=1,
	RC_BAD_EVENT=2,
	RC_DCTOOLS_ERROR=3,
	RC_STOP=-1
};

//track parameters DC
static double Rc, Xc, Yc, Zc, Za, X0, Y0, Z0, Ux, Uy, Uz, Lph;
static double Ph1, Ph2;
static double zIn, zOut, phiIn, phiOut, rIn, rOut;
//error or warning description string
static string message;
static DCTrackExt ext;

static void addintersection(int t,int cnt,int ncall,double rin,double rout,double phiin,double phiout,
	double zin,double zout,double ph1,double ph2,double dx, double dy, double dz)
{
	double lr,lphi,lz;
	int icnt=cnt-1;

	//link counter to track
	int tlast=atc_track.cnt_ntrk[icnt];
	atc_track.cnt_tracks[icnt][tlast]=t+1;
	atc_track.cnt_ntrk[icnt]++;

	//link track to counter
	int last=atc_track.ncnt_on_track[t];         //current counter index for the track
	atc_track.ncnt_on_track[t]++;
	atc_track.cnt_cross[t][last]=cnt;

	atc_track.track_end[t][last]=ncall;
        //get inlet position for track in counter
	atc_get_local(cnt,rin,phiin,zin,lr,lphi,lz);
	atc_track.rin[t][last]=lr;
	atc_track.phiin[t][last]=lphi;
	atc_track.zin[t][last]=lz;
	atc_track.ph_in[t][last]=ph1;
        //get exit position for track in counter
	atc_get_local(cnt,rout,phiout,zout,lr,lphi,lz);
	atc_track.rout[t][last]=lr;
	atc_track.phiout[t][last]=lphi;
	atc_track.zout[t][last]=lz;
	atc_track.ph_out[t][last]=ph2;

        //lenght track in aerogel
	atc_track.tlen_in_aer[t][last]=Lph*fabs(ph2-ph1);
	double lphi1=atc_track.phiin[t][last], lphi2=atc_track.phiout[t][last];

	//check if track cross WLS and get track length in WLS
	if( lphi1*lphi2 < 0 ) {
		atc_track.near_wls[t][last]=1;
		atc_track.wls_hit[t][last]=1;
		double phi=atc_geom.phi[cnt-1];
		double ct; // cosine of angle between the track direction and the WLS normal
		double cl; // cosine of angle between the track direction and WLS longitudinal direction
		double ph, r, z;
		if( ext.IntersectZplane(phi,ph1,ph2,dx,dy,dz,ph,r,z) ) { //possibly rounding error occured
			ph=(ph1+ph2)/2;
			r=(rin+rout)/2;
			z=(zin+zout)/2;
		}
		atc_track.rwls[t][last]=r;
		atc_track.phiwls[t][last]=phi;
		atc_track.zwls[t][last]=z;
		if( fMagField ) {                      // helix track
			ct=Rc*fabs(sin(ph)*sin(phi)+cos(ph)*cos(phi))/Lph;
			cl=Za/Lph;
		} else {                               // straight track
			ct=fabs(Ux*sin(phi)-Uy*cos(phi))/Lph;
			cl=Uz/Lph;
		}
		atc_track.tlen_in_wls[t][last] = (3./ct)<?(60./cl); //minimum of two numbers
	}

	//check if track goes near WLS: +-(3*sigma_dc+2) mm
	if( fabs(lphi1*rin)<3*atc_data.dcSigmaT+2 || fabs(lphi2*rout)<3*atc_data.dcSigmaT+2 )
		atc_track.near_wls[t][last]=1;
}

//track hit in ATC
static int track2atc(int t,int ncall=0)
{
	static char str[100];

	bool hitInBarrel=false, hitInEndcap=false;
	bool extendToEndcap_In11=true, left=true;
	bool extendToEndcap_In12=true;

	bool hitOutBarrel=false, hitOutEndcap=false;
	bool extendToEndcap_Out21=true;
	bool extendToEndcap_Out22=true;

	// intresection points
	int cnt1, cnt2;
	double zMin, zMax, r, z, phi, dummy;
	double phaseIn, phaseOut, ph;

	message.clear();

        int eventNumber=kedrraw_.Header.Number;

        //coordinates input/output for ATC counter
	zIn=zOut=phiIn=phiOut=rIn=rOut=0;
	phaseIn=phaseOut=0;

	if( fMagField ) {
		ext.linear=false;
		ext.x0=Xc; ext.y0=Yc; ext.z0=Zc;
		ext.rc=Rc; ext.za=Za;
		left=(Za*(Ph2-Ph1)>0);
	} else {
                //not mag. field - linear
	        ext.linear=true;
		ext.x0=X0; ext.y0=Y0; ext.z0=Z0;
		ext.vx=Ux; ext.vy=Uy; ext.vz=Uz;
		left=(Uz*(Ph2-Ph1)>0);
	}
	ext.phase1=Ph1; ext.phase2=Ph2;

        if( left ) {
            dx1=dx_el; dy1=dy_el; dz1=dz_el;
	}
	else {
            dx1=dx_er; dy1=dy_er; dz1=dz_er;
	}

//1 layer - cross inner barrel surface
if( ext.IntersectCylinder(rIn=RinBarrel1,ext.phase2,dx,dy,dz,phaseIn,zIn,phiIn)==0 )
{
    if( zIn>ZbarrelMin && zIn<ZbarrelMax ) {                                                          //actually hits the barrel
  	hitInBarrel=true;
     	if( ext.IntersectCylinder(rOut=RoutBarrel1,phaseIn,dx,dy,dz,phaseOut,zOut,phiOut)==0 ) {      //cross outer barrel surface
	    if( zOut>ZbarrelMax || zOut<ZbarrelMin ) {
					zOut=(zOut>ZbarrelMax)?ZbarrelMax:ZbarrelMin;
					ext.IntersectXYplane(zOut,dx,dy,dz,phaseOut,rOut,phiOut);
				}
			}
		}
           if( ext.IntersectCylinder(RmaxEndcap,phaseOut,dx1,dy1,dz1,ph,z,phi)==0 ) {
			if( z>ZinRightEndcap1 && z<ZinLeftEndcap1 )                           // actually goes away
			extendToEndcap_In11=false;                                            //do not continue with endcap
			if( z>ZinRightEndcap2 && z<ZinLeftEndcap2 )                           // actually goes away
			extendToEndcap_Out21=false;                                           //do not continue with endcap
		} else {
			message="track looping not implemented";
			return RC_BAD_EVENT;
		}
} //end if extend to barrel

//there is hit in innner barrel
if( hitInBarrel )
{
		int nc=ncall*2;
		atc_track.rB[t][nc]=rIn; atc_track.rB[t][nc+1]=rOut;
		atc_track.phiB[t][nc]=phiIn; atc_track.phiB[t][nc+1]=phiOut;
		atc_track.zB[t][nc]=zIn; atc_track.zB[t][nc+1]=zOut;

		cnt1=atc_get_barrel_in_counter(zIn,phiIn);
		cnt2=atc_get_barrel_in_counter(zOut,phiOut);
}
if( hitInBarrel )
{
		if( cnt1!=cnt2 ) { // two or more counters crossed
			int nb = atc_neighbours(cnt1,cnt2);             //check - are counters neighbours
			if( abs(nb)==1 ) {                              //if the counters are neighbours on phi
				phi=atc_geom.phi[cnt1-1]-HalfPhiPitch*nb;
				if( ext.IntersectZplane(phi,phaseIn,phaseOut,dx,dy,dz,ph,r,z) ) { // rounding error
					sprintf(str,"rounding problem in hitting barrel phi-adjacent counters: dPh=%5.3g",phaseOut-phaseIn);
					message=str;
					ph=(phaseIn+phaseOut)/2; r=(rIn+rOut)/2; z=(zIn+zOut)/2;
				}
			} else if( abs(nb)==2 ) {
				if( ext.IntersectXYplane(z=ZboundaryBarrel1,dx,dy,dz,ph,r,phi) ) {
					sprintf(str,"rounding problem in hitting barrel z-adjacent counters: dPh=%5.3g",phaseOut-phaseIn);
					message=str;
					ph=(phaseIn+phaseOut)/2; r=(rIn+rOut)/2; phi=(phiIn+phiOut)/2;
				}
			}
			else { // nb==0 - not neighbours
				sprintf(str,"more than 2 counters crossed in barrel: in %d, out %d",cnt1,cnt2);
				message=str;
				message+="\nNot implemented yet.";
				return RC_BAD_EVENT;
			}
                        addintersection(t,cnt1,ncall,rIn,r,phiIn,phi,zIn,z,phaseIn,ph,dx,dy,dz);
			addintersection(t,cnt2,ncall,r,rOut,phi,phiOut,z,zOut,ph,phaseOut,dx,dy,dz);
		} else { //one counter crossed
			addintersection(t,cnt1,ncall,rIn,rOut,phiIn,phiOut,zIn,zOut,phaseIn,phaseOut,dx,dy,dz);
		}
}

//external layer barrel
       	zIn=zOut=phiIn=phiOut=rIn=rOut=0;
        phaseIn=phaseOut=0;

        double X;

	if( hitInBarrel ) X=ext.phase2;
        else X=phaseOut;

if( ext.IntersectCylinder(rIn=RinBarrel2,X,dx,dy,dz,phaseIn,zIn,phiIn)==0 ) {                            //cross out barrel surface
      	    if( zIn>ZbarrelMin && zIn<ZbarrelMax ) {                                                     //actually hits the barrel
		        hitOutBarrel=true;                                                               //hit in external barrel
			if( ext.IntersectCylinder(rOut=RoutBarrel2,phaseIn,dx,dy,dz,phaseOut,zOut,phiOut)==0 ) {
			    if( zOut>ZbarrelMax || zOut<ZbarrelMin ) {                                   //but endcap earlier
				    zOut=(zOut>ZbarrelMax)?ZbarrelMax:ZbarrelMin;
					ext.IntersectXYplane(zOut,dx,dy,dz,phaseOut,rOut,phiOut);
				}
			}
		}
           if( ext.IntersectCylinder(RmaxEndcap,phaseOut,dx1,dy1,dz1,ph,z,phi)==0 ) {
			if( z>ZinRightEndcap1 && z<ZinLeftEndcap1 )                           //actually goes away
			extendToEndcap_In12=false;                                            //do not continue with endcap
			if( z>ZinRightEndcap2 && z<ZinLeftEndcap2 )                           //actually goes away
	        	extendToEndcap_Out22=false;                                           //do not continue with endcap
		}
	  else {
			message="track looping not implemented";
			return RC_BAD_EVENT;
		}

} //end if extend to barrel

//there is hit in external barrel
if( hitOutBarrel )
{
		int nc=ncall*2;
        	atc_track.rB[t][nc]=rIn; atc_track.rB[t][nc+1]=rOut;
		atc_track.phiB[t][nc]=phiIn; atc_track.phiB[t][nc+1]=phiOut;
		atc_track.zB[t][nc]=zIn; atc_track.zB[t][nc+1]=zOut;

		cnt1=atc_get_barrel_out_counter(zIn,phiIn);
		cnt2=atc_get_barrel_out_counter(zOut,phiOut);
}
if( hitOutBarrel )
{

		if( cnt1!=cnt2 ) { // two or more counters crossed
			int nb = atc_neighbours(cnt1,cnt2);
			if( abs(nb)==1 ) {
				phi=atc_geom.phi[cnt1-1]-HalfPhiPitch*nb;
				if( ext.IntersectZplane(phi,phaseIn,phaseOut,dx,dy,dz,ph,r,z) ) { // rounding error
					sprintf(str,"rounding problem in hitting barrel phi-adjacent counters: dPh=%5.3g",phaseOut-phaseIn);
					message=str;
					ph=(phaseIn+phaseOut)/2; r=(rIn+rOut)/2; z=(zIn+zOut)/2;
				}
			} else if( abs(nb)==2 ) {
				if( ext.IntersectXYplane(z=ZboundaryBarrel2,dx,dy,dz,ph,r,phi) ) {
					sprintf(str,"rounding problem in hitting barrel z-adjacent counters: dPh=%5.3g",phaseOut-phaseIn);
					message=str;
					ph=(phaseIn+phaseOut)/2; r=(rIn+rOut)/2; phi=(phiIn+phiOut)/2;
				}
			}
			else { // nb==0 - not neighbours
				sprintf(str,"more than 2 counters crossed in barrel: in %d, out %d",cnt1,cnt2);
				message=str;
				message+="\nNot implemented yet.";
				return RC_BAD_EVENT;
			}
			addintersection(t,cnt1,ncall,rIn,r,phiIn,phi,zIn,z,phaseIn,ph,dx,dy,dz);
			addintersection(t,cnt2,ncall,r,rOut,phi,phiOut,z,zOut,ph,phaseOut,dx,dy,dz);
		} else {                                                                              //one counter crossed
			addintersection(t,cnt1,ncall,rIn,rOut,phiIn,phiOut,zIn,zOut,phaseIn,phaseOut,dx,dy,dz);
		}
}
//cross inner endcap
	if( !extendToEndcap_In11 ) return 0;         //track does not unput to endcap
	if( !extendToEndcap_In12 ) return 0;

        //track hit in left endcap
	if( left ) {
	    zMin=zIn=ZinLeftEndcap1; zMax=zOut=ZoutLeftEndcap1;
            dx1=dx_el; dy1=dy_el; dz1=dz_el;
	}
	else {                                                      //track hit in right endcap
	   zMax=zIn=ZinRightEndcap1; zMin=zOut=ZoutRightEndcap1;
           dx1=dx_er; dy1=dy_er; dz1=dz_er;
	}

	ext.IntersectXYplane(zIn,dx1,dy1,dz1,phaseIn,rIn,phiIn);
	ext.IntersectXYplane(zOut,dx1,dy1,dz1,phaseOut,rOut,phiOut);

	if( rIn<RminEndcap || rOut<RminEndcap ) {
            if( ext.IntersectCylinder(RminEndcap,ext.phase2,dx1,dy1,dz1,ph,z,phi)==0 ) {
			if( z>=zMin && z<=zMax ) {
				if( rIn<RminEndcap ) {
					rIn=RminEndcap; zIn=z;  phiIn=phi; phaseIn=ph;
				} else {
					rOut=RminEndcap; zOut=z; phiOut=phi; phaseOut=ph;
				}
				hitInEndcap=true;
			}
		}
	}else {
		if( rIn>RmaxEndcap && rOut>RmaxEndcap ) { //indeed there are some events
			message="should enter the endcap";
		        return RC_BAD_EVENT;
		}
		if( rOut>RmaxEndcap ) { //at endcap edge
			ext.IntersectCylinder(RmaxEndcap,ext.phase2,dx1,dy1,dz1,ph,z,phi);
			if( z>=zMin && z<=zMax ) {
                        	rOut=RmaxEndcap; zOut=z; phiOut=phi; phaseOut=ph;
			}
	    }
		hitInEndcap=true;
	}
//there is hit in endcap
if( hitInEndcap )
{
		int nc=ncall*2;
		atc_track.rE[t][nc]=rIn; atc_track.rE[t][nc+1]=rOut;
		atc_track.phiE[t][nc]=phiIn; atc_track.phiE[t][nc+1]=phiOut;
		atc_track.zE[t][nc]=zIn; atc_track.zE[t][nc+1]=zOut;

		int cnt1=atc_get_endcap_in_counter(left,phiIn);
		int cnt2=atc_get_endcap_in_counter(left,phiOut);

		if( cnt1!=cnt2 ) {
		    int nb = atc_neighbours(cnt1,cnt2);
			if( abs(nb)==1 ) {
				phi=atc_geom.phi[cnt1-1]-HalfPhiPitch*nb;
				if( ext.IntersectZplane(phi,phaseIn,phaseOut,dx1,dy1,dz1,ph,r,z) ) {
					sprintf(str,"rounding problem in hitting adjacent endcap counters: dPh=%5.3g",phaseOut-phaseIn);
					ph=(phaseIn+phaseOut)/2; r=(rIn+rOut)/2; z=(zIn+zOut)/2;
				}
			} else { // nb==0
				sprintf(str,"more than 2 counters crossed in endcap: in %d, out %d",cnt1,cnt2);
				message=str;
				message+="\nNot implemented yet.";
				return RC_BAD_EVENT;
				//TODO
			}
			addintersection(t,cnt1,ncall,rIn,r,phiIn,phi,zIn,z,phaseIn,ph,dx1,dy1,dz1);
			addintersection(t,cnt2,ncall,r,rOut,phi,phiOut,z,zOut,ph,phaseOut,dx1,dy1,dz1);
		} else { //one counter crossed
			addintersection(t,cnt1,ncall,rIn,rOut,phiIn,phiOut,zIn,zOut,phaseIn,phaseOut,dx1,dy1,dz1);
		}
}
//cross external endcap
	if( !extendToEndcap_Out21 ) return 0;
	if( !extendToEndcap_Out22 ) return 0;

	if( left ) {
	    zMin=zIn=ZinLeftEndcap2;	    zMax=zOut=ZoutLeftEndcap2;
            dx1=dx_el; dy1=dy_el; dz1=dz_el;
	}
	else {
	    zMax=zIn=ZinRightEndcap2;      zMin=zOut=ZoutRightEndcap2;
            dx1=dx_er; dy1=dy_er; dz1=dz_er;
	}

	ext.IntersectXYplane(zIn,dx1,dy1,dz1,phaseIn,rIn,phiIn);
	ext.IntersectXYplane(zOut,dx1,dy1,dz1,phaseOut,rOut,phiOut);


	if( hitInEndcap ) X=phaseOut;
        else X=ext.phase2;

	if( rIn<RminEndcap || rOut<RminEndcap ) {
            if( ext.IntersectCylinder(RminEndcap,X,dx1,dy1,dz1,ph,z,phi)==0 ) {
		if( z>=zMin && z<=zMax ) {
				if( rIn<RminEndcap ) {
					rIn=RminEndcap; zIn=z;  phiIn=phi; phaseIn=ph;
				} else {
					rOut=RminEndcap; zOut=z; phiOut=phi; phaseOut=ph;
				}
				hitOutEndcap=true;
			}
	    }
	} else {
		if( rIn>RmaxEndcap && rOut>RmaxEndcap ) { //indeed there are some events
			message="should enter the endcap";
			return RC_BAD_EVENT;
		}
		if( rOut>RmaxEndcap ) { //at endcap edge
			ext.IntersectCylinder(RmaxEndcap,X,dx1,dy1,dz1,ph,z,phi);
			if( z>=zMin && z<=zMax ) {
                        	rOut=RmaxEndcap; zOut=z; phiOut=phi; phaseOut=ph;
			}
		}
                hitOutEndcap=true;
	}
	if( hitOutEndcap ) {
		int nc=ncall*2;
		atc_track.rE[t][nc]=rIn; atc_track.rE[t][nc+1]=rOut;
		atc_track.phiE[t][nc]=phiIn; atc_track.phiE[t][nc+1]=phiOut;
		atc_track.zE[t][nc]=zIn; atc_track.zE[t][nc+1]=zOut;

		int cnt1=atc_get_endcap_out_counter(left,phiIn);
		int cnt2=atc_get_endcap_out_counter(left,phiOut);

		if( cnt1!=cnt2 ) { //two or more counters crossed
			int nb = atc_neighbours(cnt1,cnt2);
			if( abs(nb)==1 ) {
				phi=atc_geom.phi[cnt1-1]-HalfPhiPitch*nb;
				if( ext.IntersectZplane(phi,phaseIn,phaseOut,dx1,dy1,dz1,ph,r,z) ) {
					sprintf(str,"rounding problem in hitting adjacent endcap counters: dPh=%5.3g",phaseOut-phaseIn);
					ph=(phaseIn+phaseOut)/2; r=(rIn+rOut)/2; z=(zIn+zOut)/2;
				}
			} else { // nb==0
				sprintf(str,"more than 2 counters crossed in endcap: in %d, out %d",cnt1,cnt2);
				message=str;
				message+="\nNot implemented yet.";
				return RC_BAD_EVENT;
				//TODO
			}
			addintersection(t,cnt1,ncall,rIn,r,phiIn,phi,zIn,z,phaseIn,ph,dx1,dy1,dz1);
			addintersection(t,cnt2,ncall,r,rOut,phi,phiOut,z,zOut,ph,phaseOut,dx1,dy1,dz1);
		} else { //one counter crossed
		    addintersection(t,cnt1,ncall,rIn,rOut,phiIn,phiOut,zIn,zOut,phaseIn,phaseOut,dx1,dy1,dz1);
		}
	}
	return 0;
}

static void printtrack(int t)
{
	cerr.precision(4);
	if( fMagField ) {
		cerr<<"P="<<tP(t)<<" Xc="<<Xc<<" Yc="<<Yc<<" Zc="<<Zc<<" Rc="<<Rc<<" Za="<<Za<<"\n"
			<<"Ph0="<<tPh0(t)<<" Ph1="<<Ph1<<" Ph2="<<Ph2
			<<" X1="<<10*tX1(t)<<" Y1="<<10*tY1(t)<<" Z1="<<10*tZ1(t)
			<<"; X2="<<10*tX2(t)<<" Y2="<<10*tY2(t)<<" Z2="<<10*tZ2(t)<<endl;
	} else {
		cerr<<"Ux="<<Ux<<" Uy="<<Uy<<" Uz="<<Uz<<"\n"
			<<"X0="<<X0<<" Y0="<<Y0<<" Z0="<<Z0<<endl;
	}
	cerr.precision(0);
}

static int beam_rec(long eventNumber)                                       //events of experiment
{
        //number tracks
        atc_track.ntracks = eTracksAll>ATC_NTRK?ATC_NTRK:eTracksAll;

	for(int t=0; t<atc_track.ntracks; t++ )
	{
		if( atc_droptr[t] ) continue;                               //this track is requested to skip
		//initial track parameters
		if( fMagField ) {
			Rc=10*tRc(t);
			Xc=10*tXc(t);
			Yc=10*tYc(t);
			Zc=10*tZc(t);
			Za=10*tZa(t);
			Lph=sqrt(Rc*Rc+Za*Za);
			Ph1=tPh1(t); Ph2=tPh2(t);
		} else { //no magnetic field, use alternative parametrization
			X0=t0X(t); Y0=t0Y(t); Z0=t0Z(t);
			Ux=tUx(t);
			Uy=tUy(t);
			Uz=tUz(t);
			Lph=sqrt(1+Uz*Uz);
			Ph1=t1phase(t); Ph2=t2phase(t);
		}
		int rc=track2atc(t);
		if( !rc && !message.empty() ) {
			cerr<<" WARNING: "<<message<<endl;
			printtrack(t);
		}

		if( rc ) {
			cerr<<"ATCREC: Event "<<eventNumber<<" track "<<t;
			switch( rc ) {
				case RC_DCREC_ERROR:
					cerr<<" DC ERROR: "<<message<<endl;
					printtrack(t);
					return 1; //skip the event
					break;
				case RC_BAD_EVENT:
					cerr<<" DROP EVENT: "<<message<<endl;
					return 1; //skip the event
					break;
				case RC_DCTOOLS_ERROR:
					cerr<<" track extension ERROR: "<<message<<endl;
					printtrack(t);
					atc_track.illtrack[t]=1;
					break;
				case RC_STOP:
					cerr<<" STOP: "<<message<<endl;
					return -1;
					break;
				default:
					return 1;
			}
		}
	} // end of track loop
	memset(atc_droptr,0,ATC_NTRK);

	return 0;
}

//events of cosmic rays
static int cosmic_rec(long eventNumber)
{
	atc_track.ntracks = eTracksAll>ATC_NTRK?ATC_NTRK:eTracksAll;

	for(int t=0; t<atc_track.ntracks; t++ )
	{
		if( atc_droptr[t] ) continue; //this track is requested to skip

		if( fMagField ) {
			Rc=10*tRc(t);
			Xc=10*tXc(t);
			Yc=10*tYc(t);
			Zc=10*tZc(t);
			Za=10*tZa(t);
			Lph=sqrt(Rc*Rc+Za*Za);
			Ph1=tPh2(t); Ph2=tPh1(t);
		} else {
			X0=tX0(t); Y0=tY0(t); Z0=tZ0(t);
			Ux=tUx(t);
			Uy=tUy(t);
			Uz=tUz(t);
			Lph=sqrt(1+Uz*Uz);
			Ph1=t2phase(t); Ph2=t1phase(t);
		}

		int rc=track2atc(t);
		if( !rc ) {
			if( !message.empty() ) {
				cerr<<" WARNING: "<<message<<endl;
				printtrack(t);
			}
			//reverse the track
			double phtemp=Ph1;
			Ph1=Ph2; Ph2=phtemp;

			rc=track2atc(t,1);
			if( !rc && !message.empty() ) {
				cerr<<" WARNING: "<<message<<endl;
				printtrack(t);
			}
		}

		if( rc ) {
			cerr<<"ATCREC: Event "<<eventNumber<<" track "<<t;
			switch( rc ) {
				case RC_DCREC_ERROR:
					cerr<<" DC ERROR: "<<message<<endl;
					printtrack(t);
					return 1; //skip the event
					break;
				case RC_BAD_EVENT:
					cerr<<" DROP EVENT: "<<message<<endl;
					return 1; //skip the event
					break;
				case RC_DCTOOLS_ERROR:
					cerr<<" track extension ERROR: "<<message<<endl;
					printtrack(t);
					atc_track.illtrack[t]=1;
					break;
				case RC_STOP:
					cerr<<" STOP: "<<message<<endl;
					return -1;
					break;
				default:
					return 1;
			}
		}
	} // end of track loop
	memset(atc_droptr,0,ATC_NTRK);

	return 0;
}

//reconstruction single event
int atc_event()
{
	static pair<const short*,int> dir;
	static int cosmic=-1;
	static AppEvent event;

	int eventNumber=kedrraw_.Header.Number;

	//Automatic initialization of ATC reconstruction.
	if( !atcRec ) atc_init();

	//Automatic new run determination.
	// Macro NCurRun defined in VDDCRec/kglobparam.h and updated in kdcvdrec
	if( runNumber!=NCurRun ) {
		if( atc_run(NCurRun) ) return -2;
	}

	//Clear ATC structures
	memset(&atc_track,0,sizeof(atc_track));
	memset(&atc_rec,0,sizeof(atc_rec));
	memset(&atctrackinfo,0,sizeof(atctrackinfo));

	//Wrap AtcRec as if it is called from AppFramework
	dir.first  = &kedrraw_.Word[RawFirst(SS_ATC)];          //SS_ATC - number subsystem in ReadNat (12)
	dir.second = RawLength(SS_ATC);                         //lenght write ATC in word

	event.put("AtcRawRecord",dir);                          //put raw data ATC in AtcRawRecord on event
	event.put("eventNumber",eventNumber);

	result=atcRec->event(event);

	event.clear();

	if( result.status&AppResult::ERROR )
		cerr<<"AtcRec on return: ERROR Event "<<eventNumber<<": "<<result.message<<endl;
	if( result.status&AppResult::WARNING )
		cerr<<"AtcRec on return: ERROR Event "<<eventNumber<<": "<<result.message<<endl;
	if( result.status&AppResult::LOG )
		cerr<<"AtcRec on return: Event "<<eventNumber<<": "<<result.message<<endl;

	if( result.status&AppResult::STOP ) return -1;
	if( result.status&AppResult::SKIP ) return 1;

	//Fill hit structure atc_rec
	AtcHits *hits=&atcRec->hits;
	atc_rec.eventdamage=hits->ev_type;
	atc_rec.rawDeltaT=hits->deltaT;
	atc_rec.deltaT=hits->deltaTfinal;

	AtcHits::iterator cur=hits->begin(), last=hits->end();
	AtcHit *ahit;
	int icnt, type;

	if(kedrrun_cb_.Header.RunType!=64){               //if type of run not simulation -> experiment
	for(; cur!=last; cur++) {
		icnt=cur->first-1;
		ahit=&cur->second;
		type=ahit->getType();
		atc_rec.trg[icnt]=1;
		atc_rec.htyp[icnt]=type;
		atc_rec.time[icnt]=ahit->getTfinal();
		atc_rec.amp[icnt]=ahit->getA();
		atc_rec.rawtime[icnt]=ahit->getT();

		if( type&AtcHit::Estimated ) {
			atc_rec.npe[icnt]=ahit->getNphe();
          		if( type&AtcHit::Fitted )
				atc_rec.chi2[icnt]=ahit->getChiSqr();
		}
	}
	}
        //Simulation
	else {
        atc_sim_zero_raw_hits();

	std::vector<float> npe_sim[161];
        int cnt_sim=0;
	float sum_npe[161];
        for(int i=0; i<161; i++){
	    npe_sim[i].clear();
            sum_npe[i]=0;
	}

        kscsim0_.nHits = (int)kedrmcatchits_.Long[MCATCl_Nhits-1];     //Number of hits
	int ptr = MCATC_rhl;                                           //header length in longwords

	//fill structures raw data from simulation
	for(int i=0; i<kscsim0_.nHits; i++)
	{
	    kscsim0_.cou[i] = kedrmcatchits_.Real[ptr + iatc_CNT];
	    kscsim0_.x1[i] = kedrmcatchits_.Real[ptr + ratc_X1];
	    kscsim0_.y1[i] = kedrmcatchits_.Real[ptr + ratc_Y1];
	    kscsim0_.z1[i] = kedrmcatchits_.Real[ptr + ratc_Z1];
	    kscsim0_.x2[i] = kedrmcatchits_.Real[ptr + ratc_X2];
	    kscsim0_.y2[i] = kedrmcatchits_.Real[ptr + ratc_Y2];
	    kscsim0_.z2[i] = kedrmcatchits_.Real[ptr + ratc_Z2];
	    kscsim0_.Ia_ch[i] = kedrmcatchits_.Real[ptr + ratc_Ia_ch];
	    kscsim0_.Ia_sc[i] = kedrmcatchits_.Real[ptr + ratc_Ia_sc];
	    kscsim0_.Ish_ch[i] = kedrmcatchits_.Real[ptr + ratc_Ish_ch];
	    kscsim0_.Ish_sc[i] = kedrmcatchits_.Real[ptr + ratc_Ish_sc];
	    kscsim0_.Itef_ch[i] = kedrmcatchits_.Real[ptr + ratc_Itef_ch];
	    kscsim0_.Itef_sc[i] = kedrmcatchits_.Real[ptr + ratc_Itef_sc];
	    kscsim0_.I_NAME[i] = kedrmcatchits_.Real[ptr + ratc_I_NAME];
	    kscsim0_.I_NLEVEL[i] = kedrmcatchits_.Real[ptr + ratc_I_NLEVEL];

	    cnt_sim=(int)kscsim0_.cou[i];
	    float x3=kscsim0_.Ia_ch[i];
	    float x4=kscsim0_.Itef_ch[i];
	    float x5=kscsim0_.Ish_ch[i];
	    float x6=kscsim0_.Ia_sc[i];
	    float x7=kscsim0_.Itef_sc[i];
	    float x8=kscsim0_.Ish_sc[i];
	    int x9=(int)kscsim0_.I_NAME[i];
	    int x10=(int)kscsim0_.I_NLEVEL[i];

	    float x2; //this average npe to atc_amp1 go
	    atc_amp_coef_(&cnt_sim,&x2,&x3,&x4,&x5,&x6,&x7,&x8,&x9,&x10);

	    npe_sim[cnt_sim].push_back(x2);
	    sum_npe[cnt_sim]=accumulate(npe_sim[cnt_sim].begin(), npe_sim[cnt_sim].end(), 0.0f);

	    float a3;
	    atc_amp1_(&cnt_sim,&sum_npe[cnt_sim],&a3);

	    ptr=ptr+iatc_length;  //to next hit

	    icnt=kscsim0_.cou[i]-1;
	    ahit=0;
	    type=983568;
	    atc_rec.trg[icnt]=1;
	    atc_rec.htyp[icnt]=type;
	    atc_rec.time[icnt]=0;
	    atc_rec.amp[icnt]=0;
	    atc_rec.rawtime[icnt]=0;
	    atc_rec.npe[icnt]=a3;
            atc_rec.chi2[icnt]=0;
	}

	}

	if( !atc_tracking ) return 0;

	//Determine cosmics by KDCPAR (follow it)
	if( cosmic<0 || (bool)kdcpar_.Cosmic!=(bool)cosmic ) {
		cosmic=kdcpar_.Cosmic?1:0;
		cout<<"ATCREC: Extend DC tracks with ";
		if( cosmic ) {
			event_rec=cosmic_rec;
			cout<<"cosmic_rec";
		} else {
			event_rec=beam_rec;
			cout<<"beam_rec";
		}
		cout<<endl;
	}
	//Determine helix tracks by KDCPAR (only once)
	if( fMagField<0 ) {
		if( kdcpar_.Field!=0 ) {
			fMagField=1;
			cout<<"ATCREC: Extend helix tracks"<<endl;
		} else {
			fMagField=0;
			cout<<"ATCREC: Extend straight tracks"<<endl;
		}
	}

	//Do tracks binding to counters. atc_track structure will be filled.
	int rc = event_rec(eventNumber);

	return rc;
}

void atc_sim_zero_raw_hits(){
  kscsim0_.nHits=0;
  for(int i=0; i<500; i++){
    kscsim0_.cou[i] = 0.;
    kscsim0_.x1[i] = 0.;
    kscsim0_.y1[i] = 0.;
    kscsim0_.z1[i] = 0.;
    kscsim0_.x2[i] = 0.;
    kscsim0_.y2[i] = 0.;
    kscsim0_.z2[i] = 0.;
    kscsim0_.Ia_ch[i] = 0.;
    kscsim0_.Ia_sc[i] = 0.;
    kscsim0_.Ish_ch[i] = 0.;
    kscsim0_.Ish_sc[i] = 0.;
    kscsim0_.Itef_ch[i] = 0.;
    kscsim0_.Itef_sc[i] = 0.;
    kscsim0_.I_NAME[i] = 0.;
    kscsim0_.I_NLEVEL[i] = 0.;
  }
}


