/*
 * EuroBewerter.cpp
 *
 *  Created on: May 21, 2013
 *      Author: cfdickmann
 */

#include "EuroBewerter.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "Polynom.h"
#define pi 3.141592654

namespace std {

//double AmericanOption::stammfunktion(double x){
//	double a=sigma[0]*sqrt(T);
////	double a=1;
//	return -(1.56667E-8*exp(-0.5*a*a+a*x-0.5* x*x)*pow(-9.+x,2))  /(1.86046e-6+exp(118.752/(-9.+x)));
//}
//
//double AmericanOption::stammfunktion2(double x){
//	double a=sigma[0]*sqrt(T);
////	double a=1;
//	return -(1.56667E-8*exp(-0.5*a*a+a*x-0.5* x*x)*pow(-9.+x,2))  /(1.86046e-6+exp(118.752/(-9.+x)));
//}
//

double cnd(double x)
{
  double L, K, w ;
  /* constants */
  double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
  double const a4 = -1.821255978, a5 = 1.330274429;

  L = fabs(x);
  K = 1.0 / (1.0 + 0.2316419 * L);
  w = 1.0 - 1.0 / sqrt(2 * pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

  if (x < 0 ){
    w= 1.0 - w;
  }
  return w;
}


//
//double AmericanOption::europ(double t, double T){
//	double summe=0;
//	double d_minus=(r-delta-sigma[0]*sigma[0]/2.)*T/sigma[0]/sqrt(T);
//	double d_plus=d_minus+sigma[0]*sqrt(T);
//	for(int d=0;d<D;++d)
//	{
//		double I=stammfunktion(7)-stammfunktion(sigma[0]*sqrt(T)-d_plus);
////		printf("sadads %f\n",sigma[0]*sqrt(T)-d_plus);
////		I+=stammfunktion2(0)-stammfunktion2(sigma[0]*sqrt(T)-d_plus);
//				summe+=X0[d]*exp(-delta*T)/sqrt(2*3.141592654)*I;
//	}
////	return stammfunktion(2.)-stammfunktion(1.);
//return stammfunktion(7);
//
//	double prod=1;
//	for(int d=0;d<D;++d)
//		prod*=(1-(1-pnorm(-d_minus)));
//	return summe-Strike*exp(-r*T)+Strike*exp(-t*T)*prod;//return erf(0.1);
////	return pnorm(-0.3);
//}
//

double min(double x, double y) {
	return x < y ? x : y;
}

double stammfkt(int order, double x) {
	if (order == 0)
		return sqrt(pi / 2) * erf(x / sqrt(2));
	if (order == 1)
		return -exp(-x * x / 2);
	if (order == 2)
		return -exp(-x * x / 2) * x + sqrt(pi / 2) * erf(x / sqrt(2));
	if (order == 3)
		return -exp(-x * x / 2) * (2 + x * x);
	if (order == 4)
		return -exp(-x * x / 2) * x * (3 + x * x)
				+ 3 * sqrt(pi / 2) * erf(x / sqrt(2));
	if (order == 5)
		return -exp(-x * x / 2) * (8 + 4 * pow(x, 2) + 1 * pow(x, 4));
	if (order == 6)
		return -exp(-x * x / 2) * x * (15 + 5 * x * x + x * x * x * x)
				+ 15 * sqrt(pi / 2) * erf(x / sqrt(2));

	if (order == 7)
		return -exp(-x * x / 2)
				* (48 + 24 * x * x + 6 * x * x * x * x + pow(x, 6));

	if (order == 8)
		return -exp(-x * x / 2) * x
				* (105 + 35 * x * x + 7 * x * x * x * x + pow(x, 6))
				+ 105 * sqrt(pi / 2) * erf(x / sqrt(2));

	if (order == 9)
		return -exp(-x * x / 2)
				* (384 + 192 * pow(x, 2) + 48 * pow(x, 4) + 8 * pow(x, 6)
						+ pow(x, 8));

	printf("Error 843: order falsch in stammfkt\n");
	return 0;

}

double stammfkt2(int order, double x) {
	if (order % 2 != 0) {
		double p = 0;
		double fak = 1;
		for (int nn = order - 1; nn >= 0; nn -= 2) {
			p += fak * pow(x, nn);
			fak *= nn;
		}
		return -exp(-x * x / 2) * p;
	} else {

		double p = 0;
		double fak = 1;
		for (int nn = order - 1; nn >= 0; nn -= 2) {
			p += fak * pow(x, nn);
			if (nn != 0)
				fak *= nn;
		}
		return -exp(-x * x / 2) * p + fak * sqrt(pi / 2) * erf(x / sqrt(2));;
	}
	return 0;
}

double integralExpQ(Polynom* p, double g1, double g2) {
	if (p->length == 0)
		return 0;
	double summe = 0;
	for (int order = 0; order < p->length; ++order)
		summe += p->koeff[order]* (stammfkt2(order, g2) - stammfkt2(order, g1));
	return summe;
}

double EuroBewerter::max_call(double t, double T, double* X0, int D,
		double Strike, double r, double delta, double sigma) {

//	printf("stf1: %f\n", stammfkt(0, 4.) - stammfkt(0, -2));

//	printf("stf1: %f\n",stammfkt(5,1.));
//	printf("stf2: %f\n",stammfkt2(5,1.));
//	printf("\n");
//	printf("stf1: %f\n",stammfkt(3,-1.));
//	printf("stf2: %f\n",stammfkt2(3,-1.));
//	printf("\n");
//	printf("stf1: %f\n",stammfkt(1,1.23));
//	printf("stf2: %f\n",stammfkt2(1,1.23));
//	printf("\n");
//
//	printf("stf1: %f\n",stammfkt(6,1.));
//	printf("stf2: %f\n",stammfkt2(6,1.));
//	printf("\n");
//	printf("stf1: %f\n",stammfkt(4,-1.));
//	printf("stf2: %f\n",stammfkt2(4,-1.));
//	printf("\n");
//	printf("stf1: %f\n",stammfkt(2,1.23));
//	printf("stf2: %f\n",stammfkt2(2,1.23));
//	printf("\n");
//	printf("stf1: %f\n",stammfkt(0,1.23));
//	printf("stf2: %f\n",stammfkt2(0,1.23));
//	exit(0);
T=T-t;
	double summe = 0;

	for (int d = 0; d < D; ++d) {
		double d_minus = (  log (X0[d]/Strike)+ (r - delta - sigma * sigma / 2.) * T) / (sigma * sqrt(T));
			double d_plus = d_minus + sigma * sqrt(T);
//		printf("d_minus=%f, d_plus=%f\n",d_minus,d_plus);
		Polynom ganz;
		double eins[1] = { 1. };
		ganz.set(eins, 1);
double I;
		for (int dd = 0; dd < D; ++dd)
			if (dd != d) {
//				double pp[18] = { 0.5, 0.3960985, 0, -0.061485, 0, 0.007456, 0,
//						-5.84946E-4, 0, 2.920034E-5, 0, -9.15823E-7, 0,
//						1.740319E-8, 0, -1.826093E-10, 0, 8.10495E-13 };

				double pp[10]={0.50000000000009,
				0.38567951086190133,
				0,
				-0.05010672697589501,
				0,
				0.004103597701237448,
				0,
				-1.631010612321749E-4,
				0,
				2.4428290978202304E-6
				};

				Polynom p;
//				double qq[4]={1,1,1,1};
				p.set(pp,10);
//				double II = integralExpQ(&p, -2, 4);
//				printf("\neinfach: II %f\n\n", II);

//				p.ausgeben();
//p.verschieben(-1.);
				p.verschieben(
						sigma * sqrt(T)
								+ log(X0[d] / X0[dd]) / (sigma * sqrt(T)));

				ganz.multiply_with(&p);
//printf("dplus %f\n",d_plus);
//				I = integralExpQ(&p,-d_plus, 5.);
//				I = integralExpQ(&p, -2, 4.);
//printf("Iinnen: %f\n",I);
//				ganz.ausgeben();
			}
//		printf("tehehest %f\n", ganz.auswerten(-1.));

		I = integralExpQ(&ganz, -d_plus, 5.);
//		Polynom pt;
//		double ptt[6]={0,0,0,0,0,1};
//		pt.set(ptt,6);
//		I=integralExpQ(&pt,-1,3);
//		printf("I=%f\n", I);
//		printf("assdaasd %f\n",X0[d] * exp(-delta * T) / sqrt(2 * 3.141592654) );
		summe += X0[d] * exp(-delta * T) / sqrt(2 * 3.141592654) * I;
	}
	//	return stammfunktion(2.)-stammfunktion(1.);

	double prod = 1;
	for (int d = 0; d < D; ++d){
		double d_minus = (  log (X0[d]/Strike)+ (r - delta - sigma * sigma / 2.) * T) / (sigma * sqrt(T));
		prod *= (1 - (1 - cnd(-d_minus)));
	}
//printf("prod%f\n",prod);

double zu=-Strike * exp(-r * T) + Strike * exp(-r * T) * prod;
//printf("summe%f, zu%f\n",summe,zu);
     return (summe + zu)*exp(-r*t); //return erf(0.1);
}

//double EuroBewerter::min_put(double t, double T, double* X0, int D,
//		double Strike, double r, double delta, double sigma) {
//	return max_call(double t, double T, double* X0, int D,
//			double Strike, double r, double delta, double sigma)-;
//}

EuroBewerter::EuroBewerter() {
	// TODO Auto-generated constructor stub

}

EuroBewerter::~EuroBewerter() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
