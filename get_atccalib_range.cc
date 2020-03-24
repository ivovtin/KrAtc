#include <iostream>
#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>

#include "KcSys/fortran.h"
#include "VDDCRec/kdcpar.h"
#include "VDDCRec/kglobparam.h"

extern "C" {
    int get_atccalib_range_(int &range)
    {
	if( NSimRun<=20013 ){ range=1; }
	if( NSimRun>20013 && NSimRun<=20570 ){ range=2; }
	if( NSimRun>20570 && NSimRun<=21199 ){ range=3; }
	if( NSimRun>21199 && NSimRun<=22035 ){ range=4; }
	if( NSimRun>22035 && NSimRun<=22330 ){ range=5; }
	if( NSimRun>22330 && NSimRun<=22628 ){ range=6; }
	if( NSimRun>22628 && NSimRun<=23108 ){ range=7; }
	if( NSimRun>23108 && NSimRun<=23554 ){ range=8; }
	if( NSimRun>23554 && NSimRun<=23842 ){ range=9; }
	if( NSimRun>23842 && NSimRun<=24073 ){ range=10; }
	if( NSimRun>24073 && NSimRun<=26611 ){ range=11; }
	if( NSimRun>26611 && NSimRun<=27180 ){ range=12; }
	if( NSimRun>27180 ){ range=13; }

	//cout<<"!!atc_event!!!Check ATCrec: RunNumber="<<kedrraw_.Header.RunNumber<<"!!!!!!"<<endl;
	//cout<<"\t"<<"MC ATC system -> NSimRun="<<NSimRun<<"\t"<<"ATC range calib "<<range<<endl;
    }
}
