/*
 * semi_test.cpp
 *
 *  Created on: Jan 11, 2013
 *      Author: cfdickmann
 */

#include "AmericanOption.h"
using namespace std;


AmericanOption* zeiger5;

double * semi_ergebnisse;
double ** semi_ergebnisse_grad;
double ** semi_ergebnisse_diff;
double ** semi_ergebnisse_gamma;

void* DELEGATE_semi_test(void* data) {
	zeiger5->semi_testThread(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void AmericanOption::semi_testThread(int threadnummer) {
	double erg=0;
	double diff_ergs[D];
	double gamma[D];
	double grad[D];
	double** x = DoubleFeld(N, D);
	int seed = time(NULL) +threadnummer + getpid();
	srand(seed);
	int durchlaeufe = (int)(double)(semi_testingpaths) / (double)(Threadanzahl);
	RNG generator;
	for (int m = 0; m < durchlaeufe; ++m) {
		double ** wdiff=DoubleFeld(N,D);
		for(int n=0;n<N;++n)
			for(int d=0;d<D;++d)
				wdiff[n][d]=sqrt(dt)*generator.nextGaussian();
		Pfadgenerieren(x, wdiff,0,X0);

		for (int n = 0; n < N; ++n){
			double a= semi_f(n, x[n]);
			erg += a;
			//fuer deltas
			for(int d=0;d<D;++d)
			{
				x[n][d]*=1.00001;
				diff_ergs[d]+= (semi_f(n, x[n])-a)/(0.00001*X0[d]);
				x[n][d]/=1.00001;
			}
			//fuer gammas
			for(int d=0;d<D;++d)
			{
				double s=0;
				x[n][d]*=1.00001;
				s+=semi_f(n, x[n]);
				x[n][d]/=1.00001;
				x[n][d]*=0.99999;
				s+=semi_f(n, x[n]);
				x[n][d]/=0.99999;
				s-=2.*a;
				gamma[d]+=s/pow(0.00001*X0[d],2.);
			}

			for(int d=0;d<D;++d)
				grad[d]+=semi_f_Abl(n, x[n],d)* x[n][d]/X0[d];
		}
		deleteDoubleFeld(wdiff,N,D);
	}
	semi_ergebnisse[threadnummer]=erg/(double)(durchlaeufe);
	for(int d=0;d<D;++d){
		semi_ergebnisse_diff[d][threadnummer]=diff_ergs[d]/(double)(durchlaeufe);
		semi_ergebnisse_grad[d][threadnummer]=grad[d]/(double)(durchlaeufe);
		semi_ergebnisse_gamma[d][threadnummer]=gamma[d]/(double)(durchlaeufe);
	}
	deleteDoubleFeld(x,N,D);
}

void AmericanOption::semi_testing() {
	zeiger5=this;
	printf("Testing (%d testing paths) \n",semi_testingpaths);
	semi_ergebnisse=DoubleFeld(Threadanzahl);

	semi_ergebnisse_grad=DoubleFeld(D,Threadanzahl);
	semi_ergebnisse_diff=DoubleFeld(D,Threadanzahl);
	semi_ergebnisse_gamma=DoubleFeld(D,Threadanzahl);

	int* nummern=IntFeld(Threadanzahl);
	pthread_t threads[Threadanzahl];
	for (int j = 0; j < Threadanzahl; j++){
		nummern[j]=j;
		pthread_create(&threads[j], NULL, DELEGATE_semi_test,&(nummern[j]));
	}
	for (int j = 0; j < Threadanzahl; j++)
		pthread_join(threads[j], NULL);
	deleteIntFeld(nummern,Threadanzahl);

	//if(verbose)
	for (int j = 0; j < Threadanzahl; j++)
		printf("ergebni %f\n",semi_ergebnisse[j]);

	double erg=mean(semi_ergebnisse,Threadanzahl);
	double grad[D];
	double gamma[D];
	double grad_diff[D];
	for(int d=0;d<D;++d){
		grad[d]=mean(semi_ergebnisse_grad[d],Threadanzahl);
		grad_diff[d]=mean(semi_ergebnisse_diff[d],Threadanzahl);
		gamma[d]=mean(semi_ergebnisse_gamma[d],Threadanzahl);
	}
	deleteDoubleFeld(semi_ergebnisse,Threadanzahl);
	deleteDoubleFeld(semi_ergebnisse_grad,D,Threadanzahl);
	deleteDoubleFeld(semi_ergebnisse_gamma,D,Threadanzahl);
	printf("%f\n", erg);
	printf("Deltas: ");
	for(int d=0;d<D;++d)
		printf("%f, ",grad[d]);
	printf("\n");
	printf("%f\n", erg);
	printf("Diff  : ");
	for(int d=0;d<D;++d)
		printf("%f, ",grad_diff[d]);
	printf("\n");
	printf("gamma  : ");
	for(int d=0;d<D;++d)
		printf("%.10lf, ",gamma[d]);
	printf("\n");
	ErgebnisAnhaengen(erg,(char*)"ergebnisse_semi.txt");
	ErgebnisAnhaengen(grad[0],(char*)"ergebnisse_semi_Delta1.txt");
	ErgebnisAnhaengen(grad_diff[0],(char*)"ergebnisse_semi_diffDelta1.txt");
	ErgebnisAnhaengen(gamma[0],(char*)"ergebnisse_semi_gamma11.txt");
	if(D>1){
		ErgebnisAnhaengen(grad[1],(char*)"ergebnisse_semi_Delta2.txt");
		ErgebnisAnhaengen(grad_diff[1],(char*)"ergebnisse_semi_diffDelta2.txt");
	}
}


double * koeff_ergebnisse;
double * actualkoeff;



void* DELEGATE_koeff_testen_THREAD(void* data) {
	zeiger5->koeff_testen_THREAD(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

double AmericanOption::koeff_testen(double* koeff)
{
	zeiger5=this;
	time_t time1 = time(NULL);
	actualkoeff=koeff;
	koeff_ergebnisse=DoubleFeld(Threadanzahl);
	pthread_t threads[Threadanzahl];
	int* nummern=IntFeld(Threadanzahl);
	for (int t = 0; t < Threadanzahl; t++){
		nummern[t]=t;
		pthread_create(&threads[t], NULL, DELEGATE_koeff_testen_THREAD,  &(nummern[t]));
	}
	for (int t = 0; t < Threadanzahl; t++)
		pthread_join(threads[t], NULL);
	deleteIntFeld(nummern,Threadanzahl);

	if(verbose)printf("--------------Time for koeff_testing:%ld seconds\n", time(NULL) - time1);
	double result=summe(koeff_ergebnisse,Threadanzahl);
	deleteDoubleFeld(koeff_ergebnisse,Threadanzahl);

	return result;
}



void AmericanOption::koeff_testen_THREAD(int threadnummer)
{
	double ergebnis=0;
	int* koeff_index = IntFeld(Mphi);
	int indexlauf = 0;

	if(actualkoeff!=NULL)
		for (int m = 0; m < Mphi; ++m)
			if (actualkoeff[m] != 0) {
				koeff_index[indexlauf] = m;
				indexlauf++;
			}

	for (int lauf = 0; lauf < 10000; ++lauf)
		if(lauf%Threadanzahl==threadnummer)
		{
			double sum = 0;
			if(actualkoeff!=NULL)
				for (int m = 0; m < indexlauf; ++m)
					sum += actualkoeff[koeff_index[m]] * semi_Basisfunktionen(nactual, koeff_index[m], koeff_testingpaths[lauf][nactual]);
			ergebnis +=  max(0, payoff(koeff_testingpaths[lauf][nactual], nactual) - sum)/10000.;
		}
	//printf("koeff test %d: %f\n",threadnummer,ergebnis);
	deleteIntFeld(koeff_index,Mphi);
	koeff_ergebnisse[threadnummer]=ergebnis;
}

