#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
//#include <cmath>
#include "atc_geometry.h"

using namespace std;

AtcGeometry atc_geom;

static const double k2PI=2*M_PI;

int atc_neighbours(int cnt1, int cnt2) //TODO: extend to 160 counters
{
	int tw1=(cnt1-1)/20, tw2=(cnt2-1)/20;

	if( tw1==tw2 ) {
		int diff=cnt2-cnt1;
		if( abs(diff)==1 ) return diff;
		if( abs(diff)==19 ) return diff<0?1:-1;
	} else if( tw1%4==1 && tw2%4==2 ) {
		if( cnt1+20==cnt2 ) return 2;
	} else if( tw1==2 && tw2==1 ) {
		if( cnt1==cnt2+20 ) return -2;
	}
	return 0;
}
//transformation to local coordinates of a counter (0,0,0)
int atc_get_local(int cnt,double r,double phi,double z,double& lr,double& lphi,double& lz)
{
	int i=cnt-1;
	lr=r;
	lphi=phi-atc_geom.phi[i];
	lphi+=lphi>M_PI?-k2PI:(lphi<-M_PI?k2PI:0);
	lz=z-atc_geom.z[i];
	if( cnt<=40 ) {
		lz=-lz;
		lphi=-lphi;
	}
        if( cnt>=81 && cnt<=120 ) {
		lz=-lz;
		lphi=-lphi;
	}

	return 0;
}
//transformation to global coordinates (0,0,0)
int atc_get_global(int cnt,double lr,double lphi,double lz,double& r,double& phi,double &z)
{
	int i=cnt-1;
	r=lr;
	if(cnt<=80)
     {
	phi=(cnt<=40?-lphi:lphi)+atc_geom.phi[i];
	phi+=phi<0?k2PI:(phi>k2PI?-k2PI:0);
	z=(cnt<=40?-lz:lz)+atc_geom.z[i];
     }
        if(cnt>80)
     {
        phi=(cnt<=120?-lphi:lphi)+atc_geom.phi[i];
	phi+=phi<0?k2PI:(phi>k2PI?-k2PI:0);
	z=(cnt<=120?-lz:lz)+atc_geom.z[i];
     }
	return 0;
}
//internal barrel
int atc_get_barrel_in_counter(double z,double phi)   //functions to get counter number by position
{
	int CntShift=20;
	if( z<ZboundaryBarrel1 ) CntShift=40;

	while( phi<0 || phi>k2PI ) { phi+=phi<0?k2PI:-k2PI; }

	//fraction of counter measured beginning from the edge of the first one
	double fcnt=(Phi0Barrel1-phi)/PhiPitch+0.5;
	if( fcnt<0 ) fcnt+=20;

	int cnt=int(fcnt)+1;

	if(cnt>20) {
		cerr<<"atc_get_barrel_counter DEBUG: illegal number of counter "<<cnt
			<<", fcnt="<<fcnt<<", z="<<z<<", phi/pi="<<phi/M_PI<<endl;
		cnt=20;
	}

	return CntShift+cnt;
}
//external barrel
int atc_get_barrel_out_counter(double z,double phi)   //functions to get counter number by position
{
	int CntShift=100;
	if( z<ZboundaryBarrel2 ) CntShift=120;

	while( phi<0 || phi>k2PI ) { phi+=phi<0?k2PI:-k2PI; }

	//fraction of counter measured beginning from the edge of the first one
	double fcnt=(Phi0Barrel2-phi)/PhiPitch+0.5;
	if( fcnt<0 ) fcnt+=20;

	int cnt=int(fcnt)+1;

	if(cnt>20) {
		cerr<<"atc_get_barrel_counter DEBUG: illegal number of counter "<<cnt
			<<", fcnt="<<fcnt<<", z="<<z<<", phi/pi="<<phi/M_PI<<endl;
		cnt=20;
	}

	return CntShift+cnt;
}
//internal endcap
int atc_get_endcap_in_counter(int left,double phi)
{
	int CntShift;
	if( left ) CntShift=0; else CntShift=60;

	while( phi<0 || phi>k2PI ) { phi+=phi<0?k2PI:-k2PI; }

	//fraction of counter measured beginning from the edge of the first one
	double fcnt=(Phi0Endcap1-phi)/PhiPitch+0.5;
	if( fcnt<0 ) fcnt+=20;

	int cnt=int(fcnt)+1;

	if(cnt>20) {
		cerr<<"atc_get_endcap_counter DEBUG: illegal number of counter "<<cnt
			<<", fcnt="<<fcnt<<", isleft="<<(left?"true":"false")<<", phi/pi="<<phi/M_PI<<endl;
		cnt=20;
	}

	return CntShift+cnt;
}
//external endcap
int atc_get_endcap_out_counter(int left,double phi)
{
	int CntShift;
	if( left ) CntShift=80; else CntShift=140;

	while( phi<0 || phi>k2PI ) { phi+=phi<0?k2PI:-k2PI; }

	//fraction of counter measured beginning from the edge of the first one
	double fcnt=(Phi0Endcap2-phi)/PhiPitch+0.5;
	if( fcnt<0 ) fcnt+=20;

	int cnt=int(fcnt)+1;

	if(cnt>20) {
		cerr<<"atc_get_endcap_counter DEBUG: illegal number of counter "<<cnt
			<<", fcnt="<<fcnt<<", isleft="<<(left?"true":"false")<<", phi/pi="<<phi/M_PI<<endl;
		cnt=20;
	}

	return CntShift+cnt;
}
void atc_fill_geometry()
{
	//left & right endcap  inner  -> cnt 0-20 and 61-80
	for(int i=0; i<20; i++) {
		atc_geom.r[60+i] = atc_geom.r[i] = MeanRendcap;
		atc_geom.half_r[60+i] = atc_geom.half_r[i] = Lendcap/2;
		atc_geom.phi[i] = Phi0Endcap1-i*PhiPitch;
		if( atc_geom.phi[i]<0 ) atc_geom.phi[i]+=k2PI;
		atc_geom.phi[60+i] = atc_geom.phi[i];
		atc_geom.half_phi[60+i] = atc_geom.half_phi[i] = HalfPhiSizeEndcap;
		atc_geom.z[i] = MeanZleftEndcap1;
		atc_geom.z[60+i] = MeanZrightEndcap1;
		atc_geom.half_z[60+i] = atc_geom.half_z[i] = ThicknessEndcap/2;
	}

        	//left & right endcap  out  -> cnt 81-100 and 141-160
	for(int i=80; i<100; i++) {
		atc_geom.r[60+i] = atc_geom.r[i] = MeanRendcap;
		atc_geom.half_r[60+i] = atc_geom.half_r[i] = Lendcap/2;
		atc_geom.phi[i] = Phi0Endcap2-(i-80)*PhiPitch;
		if( atc_geom.phi[i]<0 ) atc_geom.phi[i]+=k2PI;
		atc_geom.phi[60+i] = atc_geom.phi[i];
		atc_geom.half_phi[60+i] = atc_geom.half_phi[i] = HalfPhiSizeEndcap;
		atc_geom.z[i] = MeanZleftEndcap2;
		atc_geom.z[60+i] = MeanZrightEndcap2;
		atc_geom.half_z[60+i] = atc_geom.half_z[i] = ThicknessEndcap/2;
	}



	//left & right barrel  inner  -> cnt 21-60
	for(int i=20; i<40; i++) {
		atc_geom.r[20+i] = atc_geom.r[i] = MeanRbarrel1;
		atc_geom.half_r[20+i] = atc_geom.half_r[i] = ThicknessBarrel1/2;

		atc_geom.phi[i] = Phi0Barrel1-(i-20)*PhiPitch;
		if( atc_geom.phi[i]<0 ) atc_geom.phi[i]+=k2PI;
		atc_geom.phi[20+i] = atc_geom.phi[i];
		atc_geom.half_phi[20+i] = atc_geom.half_phi[i] = HalfPhiSizeBarrel1;

		atc_geom.z[i] = MeanZshortBarrel1;
		atc_geom.z[20+i] = MeanZlongBarrel1;

		atc_geom.half_z[i] = LshortBarrel/2;
		atc_geom.half_z[20+i] = LlongBarrel/2;
	}

    	//left & right barrel  out -> cnt 101-140
	for(int i=100; i<120; i++) {
		atc_geom.r[20+i] = atc_geom.r[i] = MeanRbarrel2;
		atc_geom.half_r[20+i] = atc_geom.half_r[i] = ThicknessBarrel2/2;

		atc_geom.phi[i] = Phi0Barrel2-(i-100)*PhiPitch;
		if( atc_geom.phi[i]<0 ) atc_geom.phi[i]+=k2PI;
		atc_geom.phi[20+i] = atc_geom.phi[i];
		atc_geom.half_phi[20+i] = atc_geom.half_phi[i] = HalfPhiSizeBarrel2;

                atc_geom.z[i] = MeanZlongBarrel2;
		atc_geom.z[20+i] = MeanZshortBarrel2;


		atc_geom.half_z[i] = LlongBarrel/2;
		atc_geom.half_z[20+i] = LshortBarrel/2;
	}
}

