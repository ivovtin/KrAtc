#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "DCTrackExt.h"
#include "atc_geometry.h"

using std::cout;
using std::endl;

static const double eps=1e-12;
static const double epsc=1e-2;
double epsph;

//пересечение прямой цилиндра (окружности)
int DCTrackExt::line_x_cylinder(double R,double dx, double dy,double dz,double &t,double &z,double &phi)
{
	double a,b,c,D,m,s,y,x;    //коэффициенты квадратного уравнения

      //решение квадратного уравнения  (x0+vx*t)^2 + (y0+vy*t)^2 = R^2   ищем t(время)    x^2+y^2=R^2 - уравнение окружности
      //преобразование координат трека ДК в координаты АЧС через углы Эйлера
      //x`2+y`2=R`2
      a=Ax*Ax*vx*vx+Ay*Ay*vy*vy+Az*Az*vz*vz+2*Ax*Ay*vx*vy+2*Ax*Az*vx*vz+2*Ay*Az*vy*vz;
      b=2*(Ax*Ax*vx*(x0+dx)+Ay*Ay*vy*(y0+dy)+Az*Az*vz*(z0+dz)+Ax*Ay*vx*(y0+dy)+Ax*Ay*vy*(x0+dx)+Ax*Az*vx*(z0+dz)+Ax*Az*vz*(x0+dx)+Ay*Az*vy*(z0+dz)+Ay*Az*vz*(y0+dy));
      c=Ax*Ax*(x0+dx)*(x0+dx)+Ay*Ay*(y0+dy)*(y0+dy)+Az*Az*(z0+dz)*(z0+dz)+2*Ax*Ay*(x0+dx)*(y0+dy)+2*Ax*Az*(x0+dx)*(z0+dz)+2*Ay*Az*(y0+dy)*(z0+dz)-R*R;

      if( a<eps ) return 1;
      D=b*b - 4*a*c;
      if( D<0 ) return 1;
      m=-b/(2*a);
      s=copysign(sqrt(D)/(2*a),phase2-phase1);  //copysign возвращает значение, которое имеет величину одного аргумента, а знак другого
      t=phX1=m+s;
      phX2=m-s;

      x=Ax*(x0+t*vx+dx)+Ay*(y0+t*vy+dy)+Az*(z0+t*vz+dz);
      y=Bx*(x0+t*vx+dx)+By*(y0+t*vy+dy)+Bz*(z0+t*vz+dz);
      z=Cx*(x0+t*vx+dx)+Cy*(y0+t*vy+dy)+Cz*(z0+t*vz+dz);
      phi=atan2(y,x);

      if(phi<0.) phi+=k2PI;

      return 0;
}
//пересечение спирали и цилиндра (окружности)
int DCTrackExt::helix_x_cylinder(double R,double X,double dx, double dy,double dz,double &ph,double &z,double &phi)
{
      double sPh=phase2-phase1;     //определяем направление трека - вверх-вниз

     //Ищем пересечение двух окружностей:
     //cylinder:x^2+y^2=R^2,    track:(x-x0)^2+(y-y0)^2=rc^2
     //расстояние до оси трека                от начала системы координат (0,0)
     double d=sqrt((Ax*(x0+dx)+Ay*(y0+dy)+Az*(z0+dz))*(Ax*(x0+dx)+Ay*(y0+dy)+Az*(z0+dz))+(Bx*(x0+dx)+By*(y0+dy)+Bz*(z0+dz))*(Bx*(x0+dx)+By*(y0+dy)+Bz*(z0+dz)));

     //проверим, существует ли решение
     if( d<fabs(rc-R) || d>R+rc || d<eps ) return 1;     //rc- радиус трека (спирали), R - радиус АЧС   d<fabs(rc-R) - окружность внутри другой  d>R+rc - окружности лежат отдельно

        /*       //преобразование координат
        x=Ax*(x0+rc*cos(ph)+dx)+Ay*(y0+rc*sin(ph)+dy)+Az*(z0+za*ph+dz);
        y=Bx*(x0+rc*cos(ph)+dx)+By*(y0+rc*sin(ph)+dy)+Bz*(z0+za*ph+dz);
        z=Cx*(x0+rc*cos(ph)+dx)+Cy*(y0+rc*sin(ph)+dy)+Cz*(z0+za*ph+dz);
	*/

     //K1*cosФ+K2*cos^2Ф+K3*sinФ+K4*sin^2Ф+K5*cosФ*sinФ+K6*Ф+K7*Ф^2+K8*cosФ*Ф+K9*sinФ*Ф+D=R^2
     //решение уравнения методом Ньютона
     double D,K1,K2,K3,K4,K5,K6,K7,K8,K9;
     D=(Ax*(x0+dx)+Ay*(y0+dy)+Az*(z0+dz))*(Ax*(x0+dx)+Ay*(y0+dy)+Az*(z0+dz))+(Bx*(x0+dx)+By*(y0+dy)+Bz*(z0+dz))*(Bx*(x0+dx)+By*(y0+dy)+Bz*(z0+dz));
     K1=2*Ax*Ax*(x0+dx)*rc+2*Ax*Ay*(y0+dy)*rc+2*Ax*Az*(z0+dz)*rc+2*Bx*Bx*(x0+dx)*rc+2*Bx*By*(y0+dy)*rc+2*Bx*Bz*(z0+dz)*rc;
     K2=Ax*Ax*rc*rc+Bx*Bx*rc*rc;
     K3=2*Ay*Ay*(y0+dy)*rc+2*Ax*Ay*(x0+dx)*rc+2*Ay*Az*(z0+dz)*rc+2*By*By*(y0+dy)*rc+2*Bx*By*(x0+dx)*rc+2*By*Bz*(z0+dz)*rc;
     K4=Ay*Ay*rc*rc+By*By*rc*rc;
     K5=2*Ax*Ay*rc*rc+2*Bx*By*rc*rc;
     K6=2*Az*Az*(z0+dz)*za+2*Ax*Az*(x0+dx)*za+2*Ay*Az*(y0+dy)*za+2*Bz*Bz*(z0+dz)*za+2*Bx*Bz*(x0+dx)*za+2*By*Bz*(y0+dy)*za;
     K7=Az*Az*za*za+Bz*Bz*za*za;
     K8=2*Ax*Az*rc*za+2*Bx*Bz*rc*za;
     K9=2*Ay*Az*rc*za+2*By*Bz*rc*za;

     int n=20;
     double x1[n];
     double xNext[n];

     //метод Ньютона
     for(int i=0; i<n; i++)
    {
      x1[0]=X;
      xNext[i] = x1[i]-((-R*R+D+K1*cos(x1[i])+K2*cos(x1[i])*cos(x1[i])+K3*sin(x1[i])+K4*sin(x1[i])*sin(x1[i])+K5*cos(x1[i])*sin(x1[i])+K6*x1[i]+K7*x1[i]*x1[i]+K8*cos(x1[i])*x1[i]+K9*sin(x1[i])*x1[i]))/(-K1*sin(x1[i])-2*K2*cos(x1[i])*sin(x1[i])+
      K3*cos(x1[i])+2*K4*sin(x1[i])*cos(x1[i])+K5*cos(2*x1[i])+K6+2*K7*x1[i]+K8*cos(x1[i])-K8*x1[i]*sin(x1[i])+K9*sin(x1[i])+K9*x1[i]*cos(x1[i]));

      ph=x1[i];

      epsph=((epsc/rc<epsc/za)?epsc/rc:epsc/za)>pow(10,-14)?((epsc/rc<epsc/za)?epsc/rc:epsc/za):pow(10,-14);

      if(fabs(xNext[i]-x1[i]) < epsph) break;
      x1[i+1] = xNext[i];
    }

    //азимутальный угол оси трека
    double alpha=atan2((x0+dx)*Bx+(y0+dy)*By+(z0+dz)*Bz,(x0+dx)*Ax+(y0+dy)*Ay+(z0+dz)*Az);

    //соответствующее значения фазы трека phA
    double phA=alpha;
    if( sPh>0 && phA<phase2 ) phA+=k2PI;
    else if( sPh<0 && phA>phase2 ) phA-=k2PI;

    //разность фаз соответствующей углу alpha и углу первого пересечения
    double dPhase=copysign(acos((R*R-d*d-rc*rc)/(2*d*rc)),sPh);

    //фазы пересечений
    phX1=phA-dPhase;                                   //пересечение при продолжении вперед
    phX2=phA+dPhase-copysign(k2PI,sPh);                //пересечение при продолжении назад
    //определяем продольную координату пересечений с винтовой линией
    z=Cx*(x0+rc*cos(ph)+dx)+Cy*(y0+rc*sin(ph)+dy)+Cz*(z0+za*ph+dz);

    //cylinder:x^2+y^2=R^2,    track:(x-x0)^2+(y-y0)^2=rc^2    (R*cosF-x0)^2+(R*sinF-y0)^2=rc^2  ==>  (R*cosF-Ax*(x0+dx)-Ay*(y0+dy)-Az*(z0+dz))^2+(R*sinF-Bx*(x0+dx)-By*(y0+dy)-Bz*(z0+dz))^2=rc^2

    //определяем азимутальный угол пересечений
    double dPhi=copysign(acos((R*R+d*d-rc*rc)/(2*d*R)),sPh);                     //решалось уравнение (x0-R*cos(phi))^2+(y0-R*sin(phi))^2=rc^2

    phi=alpha-dPhi; //пересечение при продолжении вперед
    if( phi<0.0 ) phi += k2PI; // берем углы от 0 до 2pi

    return 0;
}

int DCTrackExt::line_x_xyplane(double Z,double dx, double dy,double dz,double &t,double &r,double &phi)
{
	double x,y;

	if( fabs(vz)<eps ) return 1;
	//x = x0+vx*t;       //уравнение прямой
	//y = y0+vy*t;

	t = (Z-z0-dz)/vz;      //вычисляем время              //t*vz=Z-z0  ==>>  Z=z0+t*vz  ==>>   Z`=Cx*x+Cy*y+Cz*Z
	if( t*(phase2-phase1)<0 ) return 2; //нет пересечения при продолжении вперед
	phX1 = t;

	//преобразование координат трека в координаты АЧС через углы Эйлера
	x=Ax*(x0+t*vx+dx)+Ay*(y0+t*vy+dy)+Az*(z0+t*vz+dz);
	y=Bx*(x0+t*vx+dx)+By*(y0+t*vy+dy)+Bz*(z0+t*vz+dz);

	r = sqrt(x*x+y*y);
	phi = r<eps ? 0.0 : atan2(y,x);

	if( phi<0.0 ) phi += k2PI; // берем углы от 0 до 2pi

	return 0;
}

int DCTrackExt::helix_x_xyplane(double Z,double dx, double dy,double dz,double &ph,double &r,double &phi)
{
	double x,y;

	if( fabs(za)<eps ) return 1;
	//x=x0+rc*cos(ph);            //винт
	//y=y0+rc*sin(ph);

	phX1=ph=(Z-z0-dz)/za;

	//преобразование координат трека в координаты АЧС через углы Эйлера
	x=Ax*(x0+rc*cos(ph)+dx)+Ay*(y0+rc*sin(ph)+dy)+Az*(z0+za*ph+dz);
	y=Bx*(x0+rc*cos(ph)+dx)+By*(y0+rc*sin(ph)+dy)+Bz*(z0+za*ph+dz);

	if( phase1<phase2 && ph<phase1 || phase1>phase2 && ph>phase1 )  return 2; //нет пересечения вдоль трека

	r = sqrt(x*x+y*y);
	phi = r<eps ? 0.0 : atan2(y,x);

	if( phi<0.0 ) phi += k2PI;

	return 0;
}

int DCTrackExt::line_x_zplane(double angle,double dx, double dy, double dz,double &t,double& r,double& z)
{
	double s, c, x, y;

	s=sin(angle);
	c=cos(angle);

        if( fabs((Bx*vx+By*vy+Bz*vz)*c-(Ax*vx+Ay*vy+Az*vz)*s)<eps ) return 1; // нет решения

        phX1 = t = ((Ax*x0+Ax*dx+Ay*y0+Ay*dy+Az*z0+Az*dz)*s-(Bx*x0+Bx*dx+By*y0+By*dy+Bz*z0+Bz*dz)*c)/((Bx*vx+By*vy+Bz*vz)*c-(Ax*vx+Ay*vy+Az*vz)*s);   //xs=yc  -ур-е плоскости

	/*
	x = x0+vx*t;
	y = y0+vy*t;
	z = z0+vz*t;
	*/
	x=Ax*(x0+t*vx+dx)+Ay*(y0+t*vy+dy)+Az*(z0+t*vz+dz);
        y=Bx*(x0+t*vx+dx)+By*(y0+t*vy+dy)+Bz*(z0+t*vz+dz);
        z=Cx*(x0+t*vx+dx)+Cy*(y0+t*vy+dy)+Cz*(z0+t*vz+dz);

	r = sqrt(x*x+y*y);

	return 0;
}

int DCTrackExt::helix_x_zplane(double angle,double ph1,double ph2,double dx, double dy, double dz,double &ph,double& r,double& z)
{
	double phMin, phMax;
	double D,A,B,C;

	if( ph1<ph2 ) phMin=ph1, phMax=ph2;                      //selection max and min in range
	else          phMax=ph1, phMin=ph2;

        //dx=0; dy=0; dz=0;
	//решение уравнения ax=by
	//sin(angle)*x`=cos(angle)*y`  ==>> привел к виду  AcosФ+BsinФ+CФ+D=0 ==>> решаем это уравнение методом Ньютона
	D = sin(angle)*(Ax*x0+Ax*dx+Ay*y0+Ay*dy+Az*dz)-cos(angle)*(Bx*x0+Bx*dx+By*y0+By*dy+Bz*z0+Bz*dz);
	A = sin(angle)*Ax*rc-cos(angle)*Bx*rc;
	B = sin(angle)*Ay*rc-cos(angle)*By*rc;
        C = sin(angle)*Az*za-cos(angle)*Bz*za;

       int n=20;
       double x1[n];
       double kor;
       double xNext[n];
       int  count = 0;

       for(int i=0; i<n; i++)
     {
	x1[0]=phMax;
        xNext[i] = x1[i]-(-A*cos(x1[i])-B*sin(x1[i])-C*x1[i]-D)/(A*sin(x1[i])-B*cos(x1[i])-C);
        count++;

	kor=x1[i];
	epsph=((epsc/rc<epsc/za)?epsc/rc:epsc/za)>pow(10,-14)?((epsc/rc<epsc/za)?epsc/rc:epsc/za):pow(10,-14);

	if(fabs(xNext[i]-x1[i]) < epsph) break;
	x1[i+1] = xNext[i];
     }

        //проверяем попадание в отрезок [phMin,phMax]
	if( kor<=phMax )
	    ph=kor;
	else
	    return 2;

        double x,y;
	x=Ax*(x0+rc*cos(ph)+dx)+Ay*(y0+rc*sin(ph)+dy)+Az*(z0+za*ph+dz);
        y=Bx*(x0+rc*cos(ph)+dx)+By*(y0+rc*sin(ph)+dy)+Bz*(z0+za*ph+dz);
        r=sqrt(x*x+y*y);
        z=Cx*(x0+rc*cos(ph)+dx)+Cy*(y0+rc*sin(ph)+dy)+Cz*(z0+za*ph+dz);

	return 0;
}
