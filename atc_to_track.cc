#include <iostream>
#include "atc_to_track.h"
#include "atcrec.h"
#include "atc_geometry.h"
#include "atc_regions.h"

using std::cout;
using std::endl;

//struct ATC_TRACK_INFO atctrackinfo;

int atc_to_track(int t)
{
    atctrackinfo.track=t;
    atctrackinfo.ncnt=atc_track.ncnt_on_track[t];                     //число пересеченных счетчиков на трек

    int i;
    int cnt;                                                          /*list of counter numbers on the track 1..160*/

     for(i=0; i<atc_track.ncnt_on_track[t]; i++)                      //цикл по числу пересеченных счетчиков на трек
    {
      cnt=atctrackinfo.cnt[i]=atc_track.cnt_cross[t][i];              //номер счетчика который пересек трек
      atctrackinfo.npe[cnt]=atc_rec.npe[cnt-1];

      //whether track goes near or right to wls (0 or 1)
      atctrackinfo.wlshit[cnt]=atc_track.wls_hit[t][i];
      atctrackinfo.nearwls[cnt]=atc_track.near_wls[t][i];
      //geometric parameters of track crossing in mm
      atctrackinfo.tlen[cnt]=atc_track.tlen_in_aer[t][i];
      atctrackinfo.pathwls[cnt]=atc_track.tlen_in_wls[t][i];

      	if( atctrackinfo.tlen[cnt]!=0 ) {
		if( atc_is_endcap(cnt) )
			atctrackinfo.npen[cnt]=ThicknessEndcap*atctrackinfo.npe[cnt]/atctrackinfo.tlen[cnt];
		else if( atc_is_barrel_1layer(cnt) )
		    atctrackinfo.npen[cnt]=ThicknessBarrel1*atctrackinfo.npe[cnt]/atctrackinfo.tlen[cnt];
		else if( atc_is_barrel_2layer(cnt) )
		    atctrackinfo.npen[cnt]=ThicknessBarrel2*atctrackinfo.npe[cnt]/atctrackinfo.tlen[cnt];
	}

	//cout<<"ATC cnt="<<atctrackinfo.cnt[i]<<"\t"<<"Np.e.="<<atctrackinfo.npe[cnt]<<"\t"<<"Laer="<<atctrackinfo.tlen[cnt]<<
	//    "\t"<<"nNp.e.="<<atctrackinfo.npen[cnt]<<"\t"<<"wlshit="<<atctrackinfo.wlshit[cnt]<<"\t"<<"nearwls="<<atctrackinfo.nearwls[cnt]<<endl;

        //======ATC_SINGLE_CROSS - есть хотя бы одно пересечение нижней либо верхней плоскости счетчика=======================
	//если есть пересечение области аэрогеля с разным отступом (10,0,5,20 мм) от стенок и от шифтера
        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION,ATC_SINGLE_CROSS) )  atctrackinfo.single_aerogel_REGION[cnt]=1;
       	else  atctrackinfo.single_aerogel_REGION[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION0,ATC_SINGLE_CROSS) )  atctrackinfo.single_aerogel_REGION0[cnt]=1;
       	else  atctrackinfo.single_aerogel_REGION0[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION5,ATC_SINGLE_CROSS) )  atctrackinfo.single_aerogel_REGION5[cnt]=1;
       	else  atctrackinfo.single_aerogel_REGION5[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION20,ATC_SINGLE_CROSS) )  atctrackinfo.single_aerogel_REGION20[cnt]=1;
       	else  atctrackinfo.single_aerogel_REGION20[cnt]=0;

	//если есть пересечение области аэрогеля и шифтера с разным отступом (10,0,5,20 мм) от стенок
        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION,ATC_SINGLE_CROSS) )  atctrackinfo.single_active_REGION[cnt]=1;
       	else  atctrackinfo.single_active_REGION[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION0,ATC_SINGLE_CROSS) )  atctrackinfo.single_active_REGION0[cnt]=1;
       	else  atctrackinfo.single_active_REGION0[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION5,ATC_SINGLE_CROSS) )  atctrackinfo.single_active_REGION5[cnt]=1;
       	else  atctrackinfo.single_active_REGION5[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION20,ATC_SINGLE_CROSS) )  atctrackinfo.single_active_REGION20[cnt]=1;
      	else  atctrackinfo.single_active_REGION20[cnt]=0;

	if( atc_track_cross_region(t,i,ATC_CRT_REGION,ATC_SINGLE_CROSS) ) atctrackinfo.single_testreg[cnt]=1;
        else  atctrackinfo.single_testreg[cnt]=0;
	//======ATC_DOUBLE_CROSS - есть пересечение верхней и нижней плоскости счетчика(длина трека максимальна)============
	//если есть пересечение области аэрогеля с разным отступом (10,0,5,20 мм) от стенок и от шифтера
        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION,ATC_DOUBLE_CROSS) )  atctrackinfo.aerogel_REGION[cnt]=1;
       	else  atctrackinfo.aerogel_REGION[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION0,ATC_DOUBLE_CROSS) )  atctrackinfo.aerogel_REGION0[cnt]=1;
       	else  atctrackinfo.aerogel_REGION0[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION5,ATC_DOUBLE_CROSS) )  atctrackinfo.aerogel_REGION5[cnt]=1;
       	else  atctrackinfo.aerogel_REGION5[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_AEROGEL_REGION20,ATC_DOUBLE_CROSS) )  atctrackinfo.aerogel_REGION20[cnt]=1;
       	else  atctrackinfo.aerogel_REGION20[cnt]=0;

	//если есть пересечение области аэрогеля и шифтера с разным отступом (10,0,5,20 мм) от стенок
        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION,ATC_DOUBLE_CROSS) )  atctrackinfo.active_REGION[cnt]=1;
       	else  atctrackinfo.active_REGION[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION0,ATC_DOUBLE_CROSS) )  atctrackinfo.active_REGION0[cnt]=1;
       	else  atctrackinfo.active_REGION0[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION5,ATC_DOUBLE_CROSS) )  atctrackinfo.active_REGION5[cnt]=1;
       	else  atctrackinfo.active_REGION5[cnt]=0;

        if( atc_track_cross_region(t,i,ATC_ACTIVE_REGION20,ATC_DOUBLE_CROSS) )  atctrackinfo.active_REGION20[cnt]=1;
      	else  atctrackinfo.active_REGION20[cnt]=0;

	if( atc_track_cross_region(t,i,ATC_CRT_REGION,ATC_DOUBLE_CROSS) ) atctrackinfo.testreg[cnt]=1;
        else  atctrackinfo.testreg[cnt]=0;

    }

    return 0;
}
int atc_info_(int* t)
{
	return atc_to_track(*t);
}
