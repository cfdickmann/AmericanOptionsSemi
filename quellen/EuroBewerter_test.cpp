/*
 * EuroBewerter_test.cpp
 *
 *  Created on: May 21, 2013
 *      Author: cfdickmann
 */

#include "EuroBewerter.h"
#include "stdio.h"

namespace std {

static void EuroBewerterTest(){
	EuroBewerter EB;
	double X0[2]={100.,100.};
	printf("EuroBewerterTest: %.2lf (11.21) \n",EB.max_call(0,3,X0,2,100,0.05,0.1,0.2));

	double X1[3]={150.,150.};
	printf("EuroBewerterTest: %.2lf (47.18) \n",EB.max_call(0,3,X1,2,100,0.05,0.1,0.2));

	double X3[3]={70.,70.};
	printf("EuroBewerterTest: %.2lf (1.43) \n",EB.max_call(0,3,X3,2,100,0.05,0.1,0.2));

	double X4[3]={70.,110.};
	printf("EuroBewerterTest: %.2lf (9.92) \n",EB.max_call(0,3,X4,2,100,0.05,0.1,0.2));

	double X5[3]={90.,160.};
	printf("EuroBew	erterTest: %.2lf (37.25) \n",EB.max_call(0,3,X5,2,100,0.05,0.1,0.2));
}

} /* namespace std */
