#ifndef _DCTrackExt_h
#define _DCTrackExt_h

#include <cmath>

static const double k2PI=2*M_PI;

class DCTrackExt
{
public:
    bool linear; //if true track is linear else helix
	double x0, y0, z0; //origin of the linear track or center of the helix
    double vx, vy, vz; //direction of the linear track
	double rc, za; //other helix parameters
	double phase1, phase2; //phases on the track

	double phX1, phX2; //phase values got from the last intersection made

	int IntersectCylinder(double R,double X, double dx, double dy, double dz, double &ph,double &z,double &phi) {
		if( linear )
			return line_x_cylinder(R,dx,dy,dz,ph,z,phi);
		else
			return helix_x_cylinder(R,X,dx,dy,dz,ph,z,phi);
	}
	int IntersectXYplane(double Z, double dx, double dy, double dz,double &ph,double &r,double &phi) {
		if( linear )
			return line_x_xyplane(Z,dx,dy,dz,ph,r,phi);
		else
			return helix_x_xyplane(Z,dx,dy,dz,ph,r,phi);
	}
	int IntersectZplane(double angle,double ph1,double ph2,double dx, double dy, double dz,double &ph,double& r,double& z) {
		if( linear )
			return line_x_zplane(angle,dx,dy,dz,ph,r,z);
		else
			return helix_x_zplane(angle,ph1,ph2,dx,dy,dz,ph,r,z);
	}

	void GetPointDec(double ph,double& x,double& y,double& z) {
		if( linear ) {
			x=x0+vx*ph;
			y=y0+vy*ph;
			z=z0+vz*ph;
		} else {
			x=x0+rc*cos(ph);
			y=y0+rc*sin(ph);
			z=z0+za*ph;
		}
	}
	void GetPointCyl(double ph,double& r,double& phi,double& z) {
		double x, y;
		if( linear ) {
			x=x0+vx*ph;
			y=y0+vy*ph;
			z=z0+vz*ph;
		} else {
			x=x0+rc*cos(ph);
			y=y0+rc*sin(ph);
			z=z0+za*ph;
		}

		r=sqrt(x*x+y*y);
		phi = atan2(y,x);
		if( phi<0. ) phi+=k2PI;
	}

	DCTrackExt() {}
	DCTrackExt(double _x0,double _y0,double _z0,double _rc,double _za,double _ph1,double _ph2) :
		linear(false),
		x0(_x0),y0(_y0),z0(_z0),
		rc(_rc),za(_za),
		phase1(_ph1),phase2(_ph2)
	{}
	DCTrackExt(double _x0,double _y0,double _z0,double _vx,double _vy,double _vz,double _t1,double _t2) :
		linear(true),
		x0(_x0),y0(_y0),z0(_z0),
		vx(_vx),vy(_vy),vz(_vz),
		phase1(_t1),phase2(_t2)
	{}
	~DCTrackExt() {}

private:
	int line_x_cylinder(double rad,double dx, double dy, double dz,double &ph,double &z,double &phi);
	int helix_x_cylinder(double rad,double X,double dx, double dy, double dz,double &ph,double &z,double &phi);
	int line_x_xyplane(double z,double dx, double dy, double dz,double &ph,double &r,double &phi);
	int helix_x_xyplane(double z,double dx, double dy, double dz,double &ph,double &r,double &phi);
	int line_x_zplane(double angle,double dx, double dy, double dz,double &t,double& r,double& z);
	int helix_x_zplane(double angle,double ph1,double ph2,double dx, double dy, double dz,double &ph,double& r,double& z);
};

#endif
