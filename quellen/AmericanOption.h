#ifndef AMERICANOPTION_H_
#define AMERICANOPTION_H_

#include "math.h"
#include "Hilfsmittel.h"
#include "RNG.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <cstring>

#define MAX_CALL 1  //bedeutet Call im Fall D=1
#define MIN_PUT 0   //bedeutet Put im Fall D=1

#define ITO 11
#define ITOrho 12
#define EULER 13
#define CIR 14


namespace std {
class AmericanOption {
public:
	AmericanOption();
	virtual ~AmericanOption();

	int option; // MAX_CALL or MIN_PUT
	double delta; //dividend yield
	double* X0; // Spot
	double Strike; // Ausuebungspreis
	double r; // interest rate
	double* sigma; //Volatility
	double T; //Gesamtzeit
	int L;
	double eta; //for jumps
	double kappa; // for mean reversion
	double theta; // for mean reversion

	double rho;
	int Testing_Dates;
	int Training_Dates;
	int N; //time discretization
	int D;
	int M; // numer of training paths
	int Mphi;
	int K; //K=NN*5+1; //Anzahl der Basisfunktionen
	double dt;
	double*** X;

	int* Exercise_Dates;
	int number_of_Exercise_Dates;

	bool Parameter_verbose;
	bool Parameter_semi;
	bool Parameter_zehnT;
	bool Parameter_vierT;//	bool loadAlphas;
	bool Parameter_verfaelscht;
	int Threadanzahl;
	int PfadModell; //Ito or euler or CIR

	void Pfadgenerieren(double** X, int start, double * S);
	void Pfadgenerieren(double** X, double** wdiff, int start, double* S);
	void Pfadgenerieren(double** X,  int start, double* S, RNG* generator);
	void Pfadgenerieren(double** X,  int start, int ende, double* S, RNG* generator);
	void Pfadgenerieren(double** X, double** wdiff, int start,int ende, double* S);
	void Pfadgenerieren(double** X, double** wdiff);
	void Pfadgenerieren(double** X);
	void Daten();
	void neueExerciseDates(int n);
	double BoxMuller(double U1, double U2);
	void stuetzerwartung_ausrechnenThread(int k);

	double max(double d1, double d2);
	//	double nextGaussian();
	//	int Poisson(double theta);
	double payoff(double* x, int time);
	//	double* payoffAbl(double* x, int time);


//	double EuropeanPut1D_discounted(double t, double T, double S, double Strike);
//	double EuropeanCall1D_discounted(double t, double T, double S, double Strike);
//	double EuropeanPut1D(double t, double T, double S, double Strike);
//	double EuropeanCall1D(double t, double T, double S, double Strike);
//	double europeanValue(double* x, double t, double T);
//	double EuropeanOption1D_discounted(double t, double T, double S, double Strike);


	//semi members
	void  stuetzpunkte_ausrichten();
	bool* stuetzstelle_active;
	double * actualkoeff;
	double koeff_testen(double* koeff);
	void koeff_testen_THREAD(int threadnummer);
	double *** koeff_testingpaths;
	int durchlaeufe;
	double ***** semi_inner_paths;
	double** semi_betas;
	int ** semi_betas_index;
	int * semi_betas_index_max;
	double* LP_mitGLPK_Loesen(double ** Matrix, bool kleinergleich, double * RS);
	double* LP_mitGLPK_Loesen(double** Matrix, bool kleinergleich, double* RS, double* C, int Mphi,int J);

	void semi_testing();
	void semi_testThread(int threadnummer);
	int semi_testingpaths;
	void semi_mehrere_S0_testen();
	void semi_ergebnisse_ausgeben();
	void stuetzpunkte_ausgeben();
	void semi_inner_paths_erzeugen();
	void inner_paths_erzeugen_THREAD(int threadnummer);
	int durchlaufactual;

	void stuetzerwartung_ausrechnen();

	void stuetzpunkte_setzen(int n);

	void lp_ausgeben();
	int J;
	double** Matrix;
	double* C;
	double** stuetzpunkte;
	double* stuetzerwartung;
	void StuetzErwartung(int t);
	double semi_Basisfunktionen(int zeit, int j, double* x);
	//	double* semi_BasisfunktionenAbl(int zeit, int j, double* x);
	double semi_Basisfunktionen1D(int zeit, int j, double* x);
	//	double* semi_Basisfunktionen1DAbl(int zeit, int j, double* x);
	double semi_Basisfunktionen2D(int zeit, int j, double* x);
	//	double* semi_Basisfunktionen2DAbl(int zeit, int j, double* x);
	double semi_BasisfunktionenHigherD(int zeit, int j, double* x);
	//	double* semi_BasisfunktionenHigherDAbl(int zeit, int j, double* x);
	void semi();
	int nactual;
	double linearCombinationOfBasis(int zeit, double* x);
//	double* linearCombinationOfBasis_Abl(int zeit, double* x);
	double linearCombination(double* koeff, double* x);
	double semi_f(int n, double* x);
//	double semi_f_Abl(int n, double* x, int d);
};

} /* namespace std */
#endif /* AMERICANOPTION_H_ */
