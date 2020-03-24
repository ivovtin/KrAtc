#ifndef ATC_GEOMETRY_H
# define ATC_GEOMETRY_H 1

# include <math.h>
# include "KcSys/fortran.h"
# include "AtcGlobals.h"
static const double rad2deg = 180/M_PI;

///////////////////////////////////////////////////
// ATC system geometry parameters (by aerogel)
// Units:
// Dimension - mm
// Angle     - radians
//////////////////////////////////////////////////
//boundaries on aerogel

//Barrel
static const double dx=2.343545, dy=-4.0684875, dz=2.5;
static const double csi=0.00, teta=0.00, phi=0.00;
static const double Ax=(cos(csi)*cos(phi)-sin(csi)*cos(teta)*sin(phi));
static const double Ay=(-cos(csi)*sin(phi)-sin(csi)*cos(phi)*cos(phi));
static const double Az=sin(csi)*sin(teta);
static const double Bx=(sin(csi)*cos(phi)+cos(csi)*cos(teta)*sin(phi));
static const double By=(-cos(csi)*sin(teta)+cos(csi)*cos(phi)*cos(phi));
static const double Bz=-cos(csi)*sin(teta);
static const double Cx=sin(teta)*sin(phi);
static const double Cy=sin(teta)*cos(phi);
static const double Cz=cos(teta);

//Left endcap
static const double dx_el=1.426345, dy_el=0.268217, dz_el=6;
static const double csi_el=0.00, teta_el=0.00, phi_el=0.00;
static const double Ax_el=(cos(csi_el)*cos(phi_el)-sin(csi_el)*cos(teta_el)*sin(phi_el));
static const double Ay_el=(-cos(csi_el)*sin(phi_el)-sin(csi_el)*cos(phi_el)*cos(phi_el));
static const double Az_el=sin(csi_el)*sin(teta_el);
static const double Bx_el=(sin(csi_el)*cos(phi_el)+cos(csi_el)*cos(teta_el)*sin(phi_el));
static const double By_el=(-cos(csi_el)*sin(teta_el)+cos(csi_el)*cos(phi_el)*cos(phi_el));
static const double Bz_el=-cos(csi_el)*sin(teta_el);
static const double Cx_el=sin(teta_el)*sin(phi_el);
static const double Cy_el=sin(teta_el)*cos(phi_el);
static const double Cz_el=cos(teta_el);

//Right endcap
static const double dx_er=4.2864, dy_er=6.187435, dz_er=0;
static const double csi_er=0.00, teta_er=0.00, phi_er=0.00;
static const double Ax_er=(cos(csi_er)*cos(phi_er)-sin(csi_er)*cos(teta_er)*sin(phi_er));
static const double Ay_er=(-cos(csi_er)*sin(phi_er)-sin(csi_er)*cos(phi_er)*cos(phi_er));
static const double Az_er=sin(csi_er)*sin(teta_er);
static const double Bx_er=(sin(csi_er)*cos(phi_er)+cos(csi_er)*cos(teta_er)*sin(phi_er));
static const double By_er=(-cos(csi_er)*sin(teta_er)+cos(csi_er)*cos(phi_er)*cos(phi_er));
static const double Bz_er=-cos(csi_er)*sin(teta_er);
static const double Cx_er=sin(teta_er)*sin(phi_er);
static const double Cy_er=sin(teta_er)*cos(phi_er);
static const double Cz_er=cos(teta_er);

static const double ThicknessBarrel1=65, ThicknessBarrel2=69.9;
//данные по чертежам
//static const double RinBarrel1=558.5, RoutBarrel1=RinBarrel1+ThicknessBarrel1;
//static const double RinBarrel2=626.5, RoutBarrel2=RinBarrel2+ThicknessBarrel2;
//измеренные данные
static const double RinBarrel1=556.833, RoutBarrel1=RinBarrel1+ThicknessBarrel1;
static const double RoutBarrel2=697.1365, RinBarrel2=RoutBarrel2-ThicknessBarrel2;

static const double MeanRbarrel1=(RoutBarrel1+RinBarrel1)/2, MeanRbarrel2=(RoutBarrel2+RinBarrel2)/2;
static const double RinLKr=730;
static const double LbarrelTotal=1077;
static const double LshortBarrel=475.4, LlongBarrel=599.4;
static const double ZbarrelMin=-538.5, ZbarrelMax=538.5;
static const double ZboundaryBarrel1=62;
static const double ZboundaryBarrel2=-62;
static const double MeanZshortBarrel1=300.8, MeanZlongBarrel1=-238.8;
static const double MeanZshortBarrel2=-300.8, MeanZlongBarrel2=238.8;

// Endcap counters constants
static const double Lendcap=513;
static const double RminEndcap=165, RmaxEndcap=678;
static const double MeanRendcap=(RminEndcap+RmaxEndcap)/2;
static const double ThicknessEndcap=70.9;
static const double ZinLeftEndcap1=578, ZoutLeftEndcap1=ZinLeftEndcap1+ThicknessEndcap;
static const double ZinLeftEndcap2=650.4, ZoutLeftEndcap2=ZinLeftEndcap2+ThicknessEndcap;
static const double MeanZleftEndcap1=(ZinLeftEndcap1+ZoutLeftEndcap1)/2;
static const double MeanZleftEndcap2=(ZinLeftEndcap2+ZoutLeftEndcap2)/2;
static const double ZinRightEndcap1=-578, ZoutRightEndcap1=ZinRightEndcap1-ThicknessEndcap;
static const double ZinRightEndcap2=-650.4, ZoutRightEndcap2=ZinRightEndcap2-ThicknessEndcap;
static const double MeanZrightEndcap1=(ZoutRightEndcap1+ZinRightEndcap1)/2;
static const double MeanZrightEndcap2=(ZoutRightEndcap2+ZinRightEndcap2)/2;

// pitch between counters
static const double PhiPitch = M_PI/10;
static const double HalfPhiPitch = M_PI/20;

// half phi size of a counter
static const double HalfPhiSizeBarrel1 = 0.1537;
static const double HalfPhiSizeBarrel2 = 0.1543;        //((98.845-81.845)*pi/180)/2
static const double HalfPhiSizeEndcap  = 0.156;
//TODO: shifters size&position
static const double HalfPhiWLSsizeBarrel1 = 0.0027;
static const double HalfPhiWLSsizeBarrel2 = 0.0024;

// positions of the first counter in a twentys
static const double Phi0Barrel1 = M_PI/2-0.0465;
static const double Phi0Endcap1 = M_PI/2-0.035;
//static const double Phi0Endcap1 = M_PI/2-PhiPitch/2;

static const double Phi0Barrel2 = M_PI/2-PhiPitch/2-0.0465;
static const double Phi0Endcap2 = M_PI/2-PhiPitch/2-0.035;
//static const double Phi0Endcap2 = M_PI/2-PhiPitch;

struct AtcGeometry {
	float r[ATC_NCNT], phi[ATC_NCNT], z[ATC_NCNT];
	float half_r[ATC_NCNT], half_phi[ATC_NCNT], half_z[ATC_NCNT];
};

extern AtcGeometry atc_geom;

//check if the counters are neighbours (have common wall),
//returns:  0 - not neighbours,
//         +1 - phi-mates, cnt1 and cnt2 are in the clock-wise order (detector viewed from left),
//         -1 - otherwise,
//         +2 - z-mates, cnt1<cnt2,
//         -2 - otherwise.
extern int atc_neighbours(int cnt1, int cnt2);
//transformation to local coordinates of a counter
extern int atc_get_local(int cnt,double r,double phi,double z,double& lr,double& lphi,double &lz);
//transformation to global coordinates
extern int atc_get_global(int cnt,double lr,double lphi,double lz,double& r,double& phi,double &z);
//functions to get counter number by position
extern int atc_get_barrel_in_counter(double z,double phi);
extern int atc_get_barrel_out_counter(double z,double phi);
extern int atc_get_endcap_in_counter(int left,double phi); //left=1 - left endcap, left=0 - right endcap
extern int atc_get_endcap_out_counter(int left,double phi); //left=1 - left endcap, left=0 - right endcap
//function to fill atc_geom
extern void atc_fill_geometry();

inline int atc_is_barrel(int cnt) { return (cnt>20 && cnt<=60 || cnt>100 && cnt<=140)?1:0; }
inline int atc_is_short_barrel(int cnt) { return (cnt>20 && cnt<=40 || cnt>120 && cnt<=140)?1:0; }
inline int atc_is_long_barrel(int cnt) { return (cnt>40 && cnt<=60 || cnt>100 && cnt<=120)?1:0; }
inline int atc_is_left_barrel(int cnt) { return (cnt>20 && cnt<=40 || cnt>100 && cnt<=120)?1:0; }
inline int atc_is_right_barrel(int cnt) { return (cnt>40 && cnt<=60 || cnt>120 && cnt<=140)?1:0; }
inline int atc_is_endcap(int cnt) { return (cnt<=20 || cnt>60 && cnt<=100 || cnt>140)?1:0; }
inline int atc_is_left_endcap(int cnt) { return (cnt<=20 || cnt>80 && cnt<=100)?1:0; }
inline int atc_is_right_endcap(int cnt) { return (cnt>60 && cnt<=80 || cnt>140)?1:0; }
inline int atc_is_barrel_1layer(int cnt) { return (cnt>20 && cnt<=60)?1:0; }
inline int atc_is_barrel_2layer(int cnt) { return (cnt>100 && cnt<=140)?1:0; }

inline int atc_group(int cnt) { return (cnt-1)/20; }

KC_SYS_FORTRAN_NAME(atc_geom, "atc_geom");
#endif

