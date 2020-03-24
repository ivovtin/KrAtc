//нахождение попадания в определенные области счетчика
#include <iostream>
#include "atc_regions.h"
#include "AtcRegion.h"
#include "atc_geometry.h"
#include "atcrec.h"

using std::cerr;
using std::endl;
// Отступ для барельных счетчиков zInd=20,
// Отступ для торцевых счетчиков  rInd
//phiBindWLS - угол phi связанный с шифтером для барели,  phiBindWall - угол phi связанный со стенкой для барели
static double zInd=0, phiBindWLS1=10./557., phiBindWall1=10./557., phiBindWLS2=10./627., phiBindWall2=10./627.;
static double rInd=0, phiEindWLS=10./678., phiEindWall=10./678.;

static bool regions_defined=false;

//Первый слой барели и торцы
//Область аэрогеля с отступом от стенок и от шифтера
static AtcRegion AerogelE, AerogelSB1, AerogelLB1,           //объекты класса AtcRegion - торец, короткий баррельный счетчик, длинный баррельный счетчик
	AerogelE0, AerogelSB1_0, AerogelLB1_0,
	AerogelE5, AerogelSB1_5, AerogelLB1_5,
	AerogelE20, AerogelSB1_20, AerogelLB1_20;
//Активная область счетчика включая шифтер с отступом от стенок
static AtcRegion ActiveE, ActiveSB1, ActiveLB1,
    ActiveE0, ActiveSB1_0, ActiveLB1_0,
    ActiveE5, ActiveSB1_5, ActiveLB1_5,
	ActiveE20, ActiveSB1_20, ActiveLB1_20;

//Область триггера на космическом телескопе
static AtcRegion CRTtriggerE, CRTtriggerSB1, CRTtriggerLB1;


//Второй слоя барели
//Область аэрогеля с отступом от стенок и от шифтера
static AtcRegion AerogelSB2, AerogelLB2,
	AerogelSB2_0, AerogelLB2_0,
	AerogelSB2_5, AerogelLB2_5,
	AerogelSB2_20, AerogelLB2_20;
//Активная область счетчика включая шифтер с отступом от стенок
static AtcRegion ActiveSB2, ActiveLB2,
    ActiveSB2_0, ActiveLB2_0,
    ActiveSB2_5, ActiveLB2_5,
    ActiveSB2_20, ActiveLB2_20;

//Область триггера на космическом телескопе
static AtcRegion CRTtriggerSB2, CRTtriggerLB2;

static void define_atc_regions()                                  //определение областей в которых смотрю срабатывания - все в локальных координатах
{
	regions_defined=true;                                     //область в которую попала частица true
//------------------------------------
//      ТОРЦЕВОЙ СЧЕТЧИК
//------------------------------------
	phiEindWall=10./678.; phiEindWLS=10./678.;  zInd=10.;  rInd=10.;
	//Область аэрогеля в торцевом счетчике с отступом от стенок и шифтера (64%)
	// основная часть вдоль плоскости шифтера (40%)
	AerogelE.addBoxD(186+rInd,0.008+phiEindWLS,543-rInd,.1457-phiEindWall);
	// аэрогель вдоль винта и дуги до стакана ФЭУ (14%)
	AerogelE.addBoxD(543-rInd,0.0147+phiEindWLS,616-rInd,0.1540-phiEindWall);
	// "угол" около ФЭУ и электроники (10%)
	AerogelE.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);

	//Активная область торцевого счетчика с отступом от стенок (74%)
	// основная часть вдоль плоскости шифтера (64%)
	ActiveE.addBox(186+rInd,-0.1457+phiEindWall,616-rInd,0.1457-phiEindWall);
	// "угол" около ФЭУ и электроники (10%)
	ActiveE.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);

	//Области аэрогеля и активной части с разными отступами от стенок
	phiEindWall=0.;  phiEindWLS=0.;  zInd=0.;  rInd=0.;
	AerogelE0.addBoxD(186+rInd,0.008+phiEindWLS,543-rInd,.1457-phiEindWall);
	AerogelE0.addBoxD(543-rInd,0.0147+phiEindWLS,616-rInd,0.1540-phiEindWall);
	AerogelE0.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);
	ActiveE0.addBox(186+rInd,-0.1457+phiEindWall,616-rInd,0.1457-phiEindWall);
	ActiveE0.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);
	phiEindWall=5./678.; phiEindWLS=5./678.;   zInd=5.;  rInd=5.;
	AerogelE5.addBoxD(186+rInd,0.008+phiEindWLS,543-rInd,.1457-phiEindWall);
	AerogelE5.addBoxD(543-rInd,0.0147+phiEindWLS,616-rInd,0.1540-phiEindWall);
	AerogelE5.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);
	ActiveE5.addBox(186+rInd,-0.1457+phiEindWall,616-rInd,0.1457-phiEindWall);
	ActiveE5.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);
	phiEindWall=20./678.; phiEindWLS=20./678.;  zInd=20.;  rInd=20.;
	AerogelE20.addBoxD(186+rInd,0.008+phiEindWLS,543-rInd,.1457-phiEindWall);
	AerogelE20.addBoxD(543-rInd,0.0147+phiEindWLS,616-rInd,0.1540-phiEindWall);
	AerogelE20.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);
	ActiveE20.addBox(186+rInd,-0.1457+phiEindWall,616-rInd,0.1457-phiEindWall);
	ActiveE20.addBoxD(616-rInd,0.0396+phiEindWLS,678-rInd,0.1546-phiEindWall);

	//Область триггера на космическом телескопе для торцевого счетчика (9%)
	CRTtriggerE.addBoxD(366,19./416.,466,50./416.);

//--------------------------------------------
//    КОРОТКИЙ БАРЕЛЬНЫЙ СЧЕТЧИК ПЕРВОГО СЛОЯ
//--------------------------------------------
	phiBindWall1=10./557.; phiBindWLS1=10./557.;   zInd=10.;  rInd=10.;
	//Область аэрогеля в кор. барельном счетчике с отступом от стенок и шифтера (60%)
	// аэрогель около конца шифтера и разъемов (3%)
	AerogelSB1.addBoxD(-237+zInd,0.018+phiBindWLS1,-194+zInd,0.1122-phiBindWall1);
	// основная часть вдоль плоскости шифтера (50%)
	AerogelSB1.addBoxD(-194+zInd,0.003+phiBindWLS1,162-zInd,0.1537-phiBindWall1);
	// аэрогель вдоль дуги до стакана ФЭУ (4%)
	AerogelSB1.addBoxD(162-zInd,0.0152+phiBindWLS1,191-zInd,0.1537-phiBindWall1);
	// "угол" около электроники (3%)
	AerogelSB1.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);

	//Активная область в кор. барельном счетчике с отступом от стенок (70%)
	// аэрогель около конца шифтера и разъемов (3%)
	ActiveSB1.addBox(-237+zInd,-0.1122+phiBindWall1,-194+zInd,0.1122-phiBindWall1);
	// основная часть вдоль плоскости шифтера (64%)
	ActiveSB1.addBox(-194+zInd,-0.1537+phiBindWall1,191-zInd,0.1537-phiBindWall1);
	// "угол" около электроники (3%)
	ActiveSB1.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);

	//Области аэрогеля и активной части с разными отступами от стенок
	phiBindWall1=0.; phiBindWLS1=0.; zInd=0.;  rInd=0.;
	AerogelSB1_0.addBoxD(-237+zInd,0.018+phiBindWLS1,-194+zInd,0.1122-phiBindWall1);
	AerogelSB1_0.addBoxD(-194+zInd,0.003+phiBindWLS1,162-zInd,0.1537-phiBindWall1);
	AerogelSB1_0.addBoxD(162-zInd,0.0152+phiBindWLS1,191-zInd,0.1537-phiBindWall1);
	AerogelSB1_0.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);
	ActiveSB1_0.addBox(-237+zInd,-0.1122+phiBindWall1,-194+zInd,0.1122-phiBindWall1);
	ActiveSB1_0.addBox(-194+zInd,-0.1537+phiBindWall1,191-zInd,0.1537-phiBindWall1);
	ActiveSB1_0.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);
	phiBindWall1=5./557.; phiBindWLS1=5./557.; zInd=5.;  rInd=5.;
	AerogelSB1_5.addBoxD(-237+zInd,0.018+phiBindWLS1,-194+zInd,0.1122-phiBindWall1);
	AerogelSB1_5.addBoxD(-194+zInd,0.003+phiBindWLS1,162-zInd,0.1537-phiBindWall1);
	AerogelSB1_5.addBoxD(162-zInd,0.0152+phiBindWLS1,191-zInd,0.1537-phiBindWall1);
	AerogelSB1_5.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);
	ActiveSB1_5.addBox(-237+zInd,-0.1122+phiBindWall1,-194+zInd,0.1122-phiBindWall1);
	ActiveSB1_5.addBox(-194+zInd,-0.1537+phiBindWall1,191-zInd,0.1537-phiBindWall1);
	ActiveSB1_5.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);
	phiBindWall1=20./557.; phiBindWLS1=20./557.; zInd=20.;  rInd=20.;
	AerogelSB1_20.addBoxD(-237+zInd,0.018+phiBindWLS1,-194+zInd,0.1122-phiBindWall1);
	AerogelSB1_20.addBoxD(-194+zInd,0.003+phiBindWLS1,162-zInd,0.1537-phiBindWall1);
	AerogelSB1_20.addBoxD(162-zInd,0.0152+phiBindWLS1,191-zInd,0.1537-phiBindWall1);
	AerogelSB1_20.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);
	ActiveSB1_20.addBox(-237+zInd,-0.1122+phiBindWall1,-194+zInd,0.1122-phiBindWall1);
	ActiveSB1_20.addBox(-194+zInd,-0.1537+phiBindWall1,191-zInd,0.1537-phiBindWall1);
	ActiveSB1_20.addBoxD(191-zInd,0.0661+phiBindWLS1,237-zInd,0.1537-phiBindWall1);

	//Область триггера на космическом телескопе для кор. барельного счетчика (7%)
	CRTtriggerSB1.addBoxD(-83,19./591,17,50./591);

//-------------------------------------------
//    ДЛИННЫЙ БАРЕЛЬНЫЙ СЧЕТЧИК ПЕРВОГО СЛОЯ
//-------------------------------------------
	phiBindWall1=10./557.; phiBindWLS1=10./557.; zInd=10.;  rInd=10.;
	//Область аэрогеля в длин. барельном счетчике с отступом от стенок и шифтера (63%)
	// основная часть вдоль плоскости шифтера (58%)
	AerogelLB1.addBoxD(-300+zInd,0.003+phiBindWLS1,224-zInd,0.1537-phiBindWall1);
	// аэрогель вдоль дуги до стакана ФЭУ (3%)
	AerogelLB1.addBoxD(224-zInd,0.0152+phiBindWLS1,253-zInd,0.1537-phiBindWall1);
	// "угол" около электроники (3%)
	AerogelLB1.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);

	//Активная часть длин. барельного счетчика с отступом от стенок (75%)
	// основная часть вдоль плоскости шифтера (72%)
	ActiveLB1.addBox(-300+zInd,-0.1537+phiBindWall1,253-zInd,0.1537-phiBindWall1);
	// "угол" около электроники (3%)
	ActiveLB1.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);

	//Области аэрогеля и активной части с разными отступами от стенок
	phiBindWall1=0.; phiBindWLS1=0.; zInd=0.;  rInd=0.;
	AerogelLB1_0.addBoxD(-300+zInd,0.003+phiBindWLS1,224-zInd,0.1537-phiBindWall1);
	AerogelLB1_0.addBoxD(224-zInd,0.0152+phiBindWLS1,253-zInd,0.1537-phiBindWall1);
	AerogelLB1_0.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);
	ActiveLB1_0.addBox(-300+zInd,-0.1537+phiBindWall1,253-zInd,0.1537-phiBindWall1);
	ActiveLB1_0.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);
	phiBindWall1=5./557.; phiBindWLS1=5./557.;  zInd=5.;  rInd=5.;
	AerogelLB1_5.addBoxD(-300+zInd,0.003+phiBindWLS1,224-zInd,0.1537-phiBindWall1);
	AerogelLB1_5.addBoxD(224-zInd,0.0152+phiBindWLS1,253-zInd,0.1537-phiBindWall1);
	AerogelLB1_5.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);
	ActiveLB1_5.addBox(-300+zInd,-0.1537+phiBindWall1,253-zInd,0.1537-phiBindWall1);
	ActiveLB1_5.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);
	phiBindWall1=20./557.; phiBindWLS1=20./557.;  zInd=20.;  rInd=20.;
	AerogelLB1_20.addBoxD(-300+zInd,0.003+phiBindWLS1,224-zInd,0.1537-phiBindWall1);
	AerogelLB1_20.addBoxD(224-zInd,0.0152+phiBindWLS1,253-zInd,0.1537-phiBindWall1);
	AerogelLB1_20.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);
	ActiveLB1_20.addBox(-300+zInd,-0.1537+phiBindWall1,253-zInd,0.1537-phiBindWall1);
	ActiveLB1_20.addBoxD(253-zInd,0.0661+phiBindWLS1,300-zInd,0.1537-phiBindWall1);

	//Область триггера на космическом телескопе для длин. барельного счетчика (6%)
	CRTtriggerLB1.addBoxD(-145,19./591,-45,50./591);

//--------------------------------------------
//    КОРОТКИЙ БАРЕЛЬНЫЙ СЧЕТЧИК ВТОРОГО СЛОЯ
//--------------------------------------------
	phiBindWall2=10./627.; phiBindWLS2=10./627.;  zInd=10.;  rInd=10.;
	//Область аэрогеля в кор. барельном счетчике с отступом от стенок и шифтера (60%)
	// основная часть вдоль плоскости шифтера (50%)
	AerogelSB2.addBoxD(-237+zInd,0.003+phiBindWLS2,162-zInd,0.1537-phiBindWall2);
	// аэрогель вдоль дуги до стакана ФЭУ (4%)
	AerogelSB2.addBoxD(162-zInd,0.0152+phiBindWLS2,191-zInd,0.1537-phiBindWall2);
	// "угол" около электроники (3%)
	AerogelSB2.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);

	//Активная область в кор. барельном счетчике с отступом от стенок (70%)
	// основная часть вдоль плоскости шифтера (64%)
	ActiveSB2.addBox(-237+zInd,-0.1537+phiBindWall2,191-zInd,0.1537-phiBindWall2);
	// "угол" около электроники (3%)
	ActiveSB2.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);

	//Области аэрогеля и активной части с разными отступами от стенок
	phiBindWall2=0.; phiBindWLS2=0.; zInd=0.;  rInd=0.;
	//AerogelSB2_0.addBoxD(-237+zInd,0.018+phiBindWLS2,-194+zInd,0.1122-phiBindWall2);
	AerogelSB2_0.addBoxD(-237+zInd,0.003+phiBindWLS2,162-zInd,0.1537-phiBindWall2);
	AerogelSB2_0.addBoxD(162-zInd,0.0152+phiBindWLS2,191-zInd,0.1537-phiBindWall2);
	AerogelSB2_0.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);
	//ActiveSB2_0.addBox(-237+zInd,-0.1122+phiBindWall2,-194+zInd,0.1122-phiBindWall2);
	ActiveSB2_0.addBox(-237+zInd,-0.1537+phiBindWall2,191-zInd,0.1537-phiBindWall2);
	ActiveSB2_0.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);
	phiBindWall2=5./627.; phiBindWLS2=5./627.; zInd=5.;  rInd=5.;
	//AerogelSB2_5.addBoxD(-237+zInd,0.018+phiBindWLS2,-194+zInd,0.1122-phiBindWall2);
	AerogelSB2_5.addBoxD(-237+zInd,0.003+phiBindWLS2,162-zInd,0.1537-phiBindWall2);
	AerogelSB2_5.addBoxD(162-zInd,0.0152+phiBindWLS2,191-zInd,0.1537-phiBindWall2);
	AerogelSB2_5.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);
	//ActiveSB2_5.addBox(-237+zInd,-0.1122+phiBindWall2,-194+zInd,0.1122-phiBindWall2);
	ActiveSB2_5.addBox(-237+zInd,-0.1537+phiBindWall2,191-zInd,0.1537-phiBindWall2);
	ActiveSB2_5.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);
	phiBindWall2=20./627.; phiBindWLS2=20./627.; zInd=20.;  rInd=20.;
	//AerogelSB2_20.addBoxD(-237+zInd,0.018+phiBindWLS2,-194+zInd,0.1122-phiBindWall2);
	AerogelSB2_20.addBoxD(-237+zInd,0.003+phiBindWLS2,162-zInd,0.1537-phiBindWall2);
	AerogelSB2_20.addBoxD(162-zInd,0.0152+phiBindWLS2,191-zInd,0.1537-phiBindWall2);
	AerogelSB2_20.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);
	//ActiveSB2_20.addBox(-237+zInd,-0.1122+phiBindWall2,-194+zInd,0.1122-phiBindWall2);
	ActiveSB2_20.addBox(-237+zInd,-0.1537+phiBindWall2,191-zInd,0.1537-phiBindWall2);
	ActiveSB2_20.addBoxD(191-zInd,0.0661+phiBindWLS2,237-zInd,0.1537-phiBindWall2);

	//Область триггера на космическом телескопе для кор. барельного счетчика (7%)
	CRTtriggerSB2.addBoxD(-83,19./591,17,50./591);

//-------------------------------------------
//    ДЛИННЫЙ БАРЕЛЬНЫЙ СЧЕТЧИК  ВТОРОГО СЛОЯ
//-------------------------------------------
	phiBindWall2=10./627.; phiBindWLS2=10./627.; zInd=10.;  rInd=10.;
	//Область аэрогеля в длин. барельном счетчике с отступом от стенок и шифтера (63%)
	// аэрогель около конца шифтера и разъемов (3%)
        AerogelLB2.addBoxD(-300+zInd,0.018+phiBindWLS2,-257+zInd,0.1122-phiBindWall2);     //во 2 слое разъемы в длинном счетчике
	// основная часть вдоль плоскости шифтера (58%)
	AerogelLB2.addBoxD(-257+zInd,0.003+phiBindWLS2,224-zInd,0.1537-phiBindWall2);
	// аэрогель вдоль дуги до стакана ФЭУ (3%)
	AerogelLB2.addBoxD(224-zInd,0.0152+phiBindWLS2,253-zInd,0.1537-phiBindWall2);
	// "угол" около электроники (3%)
	AerogelLB2.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);

	//Активная часть длин. барельного счетчика с отступом от стенок (75%)
	// аэрогель около конца шифтера и разъемов (3%)
	ActiveLB2.addBox(-300+zInd,-0.1122+phiBindWall2,-257+zInd,0.1122-phiBindWall2);
	// основная часть вдоль плоскости шифтера (72%)
	ActiveLB2.addBox(-257+zInd,-0.1537+phiBindWall2,253-zInd,0.1537-phiBindWall2);
	// "угол" около электроники (3%)
	ActiveLB2.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);

	//Области аэрогеля и активной части с разными отступами от стенок
	phiBindWall2=0.; phiBindWLS2=0.; zInd=0.;  rInd=0.;
        AerogelLB2_0.addBoxD(-300+zInd,0.018+phiBindWLS2,-257+zInd,0.1122-phiBindWall2);
	AerogelLB2_0.addBoxD(-257+zInd,0.003+phiBindWLS2,224-zInd,0.1537-phiBindWall2);
	AerogelLB2_0.addBoxD(224-zInd,0.0152+phiBindWLS2,253-zInd,0.1537-phiBindWall2);
	AerogelLB2_0.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);
	ActiveLB2_0.addBox(-300+zInd,-0.1122+phiBindWall2,-257+zInd,0.1122-phiBindWall2);
	ActiveLB2_0.addBox(-257+zInd,-0.1537+phiBindWall2,253-zInd,0.1537-phiBindWall2);
	ActiveLB2_0.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);
	phiBindWall2=5./627.; phiBindWLS2=5./627.; zInd=5.;  rInd=5.;
        AerogelLB2_5.addBoxD(-300+zInd,0.018+phiBindWLS2,-257+zInd,0.1122-phiBindWall2);
	AerogelLB2_5.addBoxD(-257+zInd,0.003+phiBindWLS2,224-zInd,0.1537-phiBindWall2);
	AerogelLB2_5.addBoxD(224-zInd,0.0152+phiBindWLS2,253-zInd,0.1537-phiBindWall2);
	AerogelLB2_5.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);
	ActiveLB2_5.addBox(-300+zInd,-0.1122+phiBindWall2,-257+zInd,0.1122-phiBindWall2);
	ActiveLB2_5.addBox(-257+zInd,-0.1537+phiBindWall2,253-zInd,0.1537-phiBindWall2);
	ActiveLB2_5.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);
	phiBindWall2=20./627.; phiBindWLS2=20./627.; zInd=20.;  rInd=20.;
        AerogelLB2_20.addBoxD(-300+zInd,0.018+phiBindWLS2,-257+zInd,0.1122-phiBindWall2);
	AerogelLB2_20.addBoxD(-257+zInd,0.003+phiBindWLS2,224-zInd,0.1537-phiBindWall2);
	AerogelLB2_20.addBoxD(224-zInd,0.0152+phiBindWLS2,253-zInd,0.1537-phiBindWall2);
	AerogelLB2_20.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);
	ActiveLB2_20.addBox(-300+zInd,-0.1122+phiBindWall2,-257+zInd,0.1122-phiBindWall2);
	ActiveLB2_20.addBox(-257+zInd,-0.1537+phiBindWall2,253-zInd,0.1537-phiBindWall2);
	ActiveLB2_20.addBoxD(253-zInd,0.0661+phiBindWLS2,300-zInd,0.1537-phiBindWall2);

	//Область триггера на космическом телескопе для длин. барельного счетчика (6%)
	CRTtriggerLB2.addBoxD(-145,19./591,-45,50./591);
}

int atc_track_cross_region(int t, int iatc, int regtype, int crosstype)             //находим область пересеченную треком
{
	if( !regions_defined ) define_atc_regions();             //если regions_defined==true то определяем все области

	if( t<0 || t>=atc_track.ntracks ) {                                                    //если нет трека или плохой
		cerr<<__func__<<": Illegal track number for this event "<<t<<endl;
		return 0;
	}
	if( iatc<0 || iatc>=atc_track.ncnt_on_track[t] ) {                                    //если нет счетчиков на трек
		cerr<<__func__<<": Illegal ATC counter index for this event "<<iatc<<endl;
		return 0;
	}

	int cnt=atc_track.cnt_cross[t][iatc];                   //номер счетчика который пересек трек

	AtcRegion *region=0;
	if( regtype==ATC_AEROGEL_REGION ) {
		if( atc_is_endcap(cnt) )                    //==>>atc_geometry.h
			region=&AerogelE;
		else if(atc_is_barrel_1layer(cnt))
		{
		    if( atc_is_short_barrel(cnt) )
			region=&AerogelSB1;
		    else if( atc_is_long_barrel(cnt) )
			region=&AerogelLB1;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		    if( atc_is_short_barrel(cnt) )
			region=&AerogelSB2;
		    else if( atc_is_long_barrel(cnt) )
			region=&AerogelLB2;
		}
	} else if( regtype==ATC_AEROGEL_REGION0 ) {
		if( atc_is_endcap(cnt) )
			region=&AerogelE0;
		else if(atc_is_barrel_1layer(cnt))
		{
		    if( atc_is_short_barrel(cnt) )
			region=&AerogelSB1_0;
		    else if( atc_is_long_barrel(cnt) )
			region=&AerogelLB1_0;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		    if( atc_is_short_barrel(cnt) )
			region=&AerogelSB2_0;
		    else if( atc_is_long_barrel(cnt) )
			region=&AerogelLB2_0;
		}
	} else if( regtype==ATC_AEROGEL_REGION5 ) {
		if( atc_is_endcap(cnt) )
			region=&AerogelE5;
		else if(atc_is_barrel_1layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&AerogelSB1_5;
		else if( atc_is_long_barrel(cnt) )
		        region=&AerogelLB1_5;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&AerogelSB2_5;
		else if( atc_is_long_barrel(cnt) )
		        region=&AerogelLB2_5;
		}
	} else if( regtype==ATC_AEROGEL_REGION20 ) {
		if( atc_is_endcap(cnt) )
			region=&AerogelE20;
		else if(atc_is_barrel_1layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&AerogelSB1_20;
		else if( atc_is_long_barrel(cnt) )
		        region=&AerogelLB1_20;
		}
                else if(atc_is_barrel_2layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&AerogelSB2_20;
		else if( atc_is_long_barrel(cnt) )
		        region=&AerogelLB2_20;
		}
	} else if( regtype==ATC_ACTIVE_REGION ) {
		if( atc_is_endcap(cnt) )
			region=&ActiveE;
		else if(atc_is_barrel_1layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB1;
		else if( atc_is_long_barrel(cnt) )
		        region=&ActiveLB1;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB2;
		else if( atc_is_long_barrel(cnt) )
		        region=&ActiveLB2;
		}
	} else if( regtype==ATC_ACTIVE_REGION0 ) {
		if( atc_is_endcap(cnt) )
			region=&ActiveE0;
		else if(atc_is_barrel_1layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB1_0;
		else if( atc_is_long_barrel(cnt) )
		    region=&ActiveLB1_0;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB2_0;
		else if( atc_is_long_barrel(cnt) )
		        region=&ActiveLB2_0;
		}
	} else if( regtype==ATC_ACTIVE_REGION5 ) {
		if( atc_is_endcap(cnt) )
			region=&ActiveE5;
		else if(atc_is_barrel_1layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB1_5;
		else if( atc_is_long_barrel(cnt) )
		        region=&ActiveLB1_5;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB2_5;
		else if( atc_is_long_barrel(cnt) )
		        region=&ActiveLB2_5;
		}
	} else if( regtype==ATC_ACTIVE_REGION20 ) {
		if( atc_is_endcap(cnt) )
			region=&ActiveE20;
		else if(atc_is_barrel_1layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB1_20;
		else if( atc_is_long_barrel(cnt) )
        	        region=&ActiveLB1_20;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&ActiveSB2_20;
		else if( atc_is_long_barrel(cnt) )
		        region=&ActiveLB2_20;
		}
	} else if( regtype==ATC_CRT_REGION ) {
		if( atc_is_endcap(cnt) )
			region=&CRTtriggerE;
		else if(atc_is_barrel_1layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&CRTtriggerSB1;
		else if( atc_is_long_barrel(cnt) )
		        region=&CRTtriggerLB1;
		}
		else if(atc_is_barrel_2layer(cnt))
		{
		if( atc_is_short_barrel(cnt) )
			region=&CRTtriggerSB2;
		else if( atc_is_long_barrel(cnt) )
		        region=&CRTtriggerLB2;
		}
	} else {
		cerr<<__func__<<": Illegal ATC region type "<<regtype<<endl;
		return 0;
	}
	assert( region!=0 );

    float phi1=atc_track.phiin[t][iatc], phi2=atc_track.phiout[t][iatc];   //присвоение локальных координат входа/выхода в счетчик
	float x1, x2;
	if( atc_is_endcap(cnt) ) {
		x1=atc_track.rin[t][iatc];
		x2=atc_track.rout[t][iatc];
	} else {
		x1=atc_track.zin[t][iatc];
		x2=atc_track.zout[t][iatc];
	}

	if( crosstype==ATC_DOUBLE_CROSS )
		return region->isIn(x1,phi1,x2,phi2)?1:0;
	if( crosstype==ATC_SINGLE_CROSS )
		return (region->isIn(x1,phi1) || region->isIn(x2,phi2))?1:0;
	if( crosstype==ATC_IN_CROSS )
		return region->isIn(x1,phi1)?1:0;
	if( crosstype==ATC_OUT_CROSS )
		return region->isIn(x2,phi2)?1:0;

	cerr<<__func__<<": Illegal ATC crossing type "<<crosstype<<endl;
	return 0;
}

int atc_trreg_(int* nt, int* natc, int *regtype, int *crosstype)
{
	return atc_track_cross_region(*nt-1,*natc-1,*regtype,*crosstype);
}

