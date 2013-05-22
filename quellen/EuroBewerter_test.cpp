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
	double X0[3]={100.,100.,100.};
	for(int k=0;k<10000;++k)
	printf("\nEuroBewerterTest: %f\n",EB.max_call(0,3,X0,2,100,0.05,0.1,0.2));

//	printf("\nEuroBewerterTest: %f\n",EB.max_call(0,3,X0,3,100,0.05,0.1,0.2));
//	printf("\nEuroBewerterTest: %f\n",EB.max_call(0,3,X0,3,100,0.05,0.1,0.2));
//	printf("\nEuroBewerterTest: %f\n",EB.max_call(0,3,X0,3,100,0.05,0.1,0.2));
  printf("Es muss rauskommen: %f\n\n",11.21);
}


} /* namespace std */
