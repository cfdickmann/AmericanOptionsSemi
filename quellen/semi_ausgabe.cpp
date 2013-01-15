/*
 * semi_ausgabe.cpp
 *
 *  Created on: Jan 11, 2013
 *      Author: cfdickmann
 */
#include "AmericanOption.h"

using namespace std;

//double AmericanOption::semi_f_Abl(int zeit, double* x, int d) {
//
//	if (zeit == N - 1){
//		double* grad;
//		grad=payoffAbl(x, zeit);
//		double erg=grad[d];
//		deleteDoubleFeld(grad,D);
//		return erg;
//	}
//
//	double sum = 0;
//	int m;
//	for (int i = 0; i < semi_betas_index_max[zeit]; ++i) {
//		m = semi_betas_index[zeit][i];
//		sum += semi_betas[zeit][m] * semi_Basisfunktionen(zeit, m, x);
//	}
//
//	if(payoff(x, zeit) < sum)
//		return 0;
//	else
//	{
//		double* A=payoffAbl(x, zeit);
//		double sum = A[d];
//		deleteDoubleFeld(A,D);
//		int m;
//		for (int i = 0; i < semi_betas_index_max[zeit]; ++i) {
//			m = semi_betas_index[zeit][i];
//			double* B=semi_BasisfunktionenAbl(zeit, m, x);
//			sum -= semi_betas[zeit][m] * B[d];
//			deleteDoubleFeld(B,D);
//		}
//		return sum;
//	}
//}

void AmericanOption::stuetzpunkte_ausgeben()
{
	for(int j=0;j<J;++j){
		printf("stuetzpunkt:(");
		for(int d=0;d<D;++d)
			printf("%.2lf,",stuetzpunkte[j][d]);
		printf(") E: %.2lf ,payoff: %.2lf\n",stuetzerwartung[j],payoff(stuetzpunkte[j],nactual));
	}
}

void AmericanOption::semi_ergebnisse_ausgeben(){
	printf("semi_betas ausgeben \n");
	for (int j = 0; j < Mphi; ++j){
		if((Mphi-j)%1000==0)printf("\n\n");
		printf("%.0lf, ", 1000*semi_betas[nactual][j]);
	}
	printf("\n");

	//    printf("C_upperbound ausgeben \n");
	//    for (int j = 0; j < J; ++j)
	//        printf("%.3lf, ", linearCombinationOfBasis(nactual, stuetzpunkte[j], 0));
	//    printf("\n");
	//    printf("f_estimated ausgeben \n");
	//    for (int j = 0; j < J; ++j) {
	//        printf("%.3lf, ", semi_f(nactual, stuetzpunkte[j]));
	//    }
	//    printf("\n");

	double start = stuetzpunkte[0][0];
	double stop = stuetzpunkte[J - 1][0];
	{
		fstream f;
		f.open("func.data", ios::out);
		for (double t = start; t <= stop; t += 0.01) {
			double point[2];
			point[0] = t;
			f << linearCombinationOfBasis(nactual, point) << endl;
		}
		f.close();
	}

	int WJ = (int) (sqrt(J));
	fstream file_diff;
	fstream file_stue;
	fstream file_linComb;
	fstream file_Q;
	fstream file_payoff;
	fstream file_europ;
	fstream file_europnah;
	fstream file_LSM_C_estimated;
	file_payoff.open("file_payoff.data", ios::out);
	file_Q.open("file_Q.data", ios::out);
	file_diff.open("file_diff.data", ios::out);
	file_linComb.open("file_linComb.data", ios::out);
	file_stue.open("file_stue.data", ios::out);
	file_europ.open("file_europ.data", ios::out);
	file_europnah.open("file_europnah.data", ios::out);
	file_LSM_C_estimated.open("file_LSM_C_estimated.data", ios::out);

	double array_payoff[J];
	double array_Q[J];
	double array_diff[J];
	double array_linComb[J];
	double array_stue[J];
	double array_europ[J];
	double array_europnah[J];
	//double array_LSM_C_estimated[J];

	for (int j = 0; j < J; j++) {
		array_Q[j]= max(payoff(stuetzpunkte[j], nactual) -  linearCombinationOfBasis(nactual, stuetzpunkte[j]), 0);
		array_payoff[j] =payoff(stuetzpunkte[j], nactual);
		array_diff[j] =linearCombinationOfBasis(nactual, stuetzpunkte[j]) - stuetzerwartung[j] ;
		array_stue[j] =stuetzerwartung[j];
		array_linComb[j] = linearCombinationOfBasis(nactual, stuetzpunkte[j]);
		// array_europ[j] = europeanValue(stuetzpunkte[j], nactual*dt, T) ;
		// array_europnah[j] = europeanValue(stuetzpunkte[j], nactual*dt, nactual * dt + dt) ;
		//     if ((Mphi == 6 && D == 1) || (Mphi == 7 && D == 2))
		//        file_LSM_C_estimated << LSM_C_estimated(stuetzpunkte[j], nactual) << endl;
	}

	for (int j = 0; j < J; j++) {
		if (D == 2 && j % WJ == 0) {
			file_LSM_C_estimated << endl;
			file_diff << endl;
			file_linComb << endl;
			file_Q << endl;
			file_stue << endl;
			file_europ << endl;
			file_europnah << endl;
			file_payoff << endl;
		}
		file_Q <<array_Q[j] << endl;
		file_payoff <<array_payoff[j]<< endl;
		file_diff << array_diff[j] << endl;
		file_stue <<array_stue[j] << endl;
		file_linComb << array_linComb[j]<< endl;
		file_europ << array_europ[j]<< endl;
		file_europnah << array_europnah[j]<< endl;
		//     if ((Mphi == 6 && D == 1) || (Mphi == 7 && D == 2))
		//        file_LSM_C_estimated << LSM_C_estimated(stuetzpunkte[j], nactual) << endl;
	}
	file_payoff.close();
	file_diff.close();
	file_Q.close();
	file_stue.close();
	file_linComb.close();
	file_europ.close();
	file_europnah.close();
	file_LSM_C_estimated.close();
}

void AmericanOption::semi_mehrere_S0_testen() {
	printf("Error 40");
	/*
  //Europ schreiben
    {
        fstream f;
        f.open("europ.data", ios::out);
        for (double start = 60; start < 100; start += 5)

            f << EuropeanPut1D_discounted(0, T, start, Strike) << endl;

        f.close();
    }

    //Value schreiben
    fstream f2;
    f2.open("values.data", ios::out);

    for (double start = 60; start < 100; start += 5) {
        int ergebnispipe[Threadanzahl][2];
        for (int z = 0; z < Threadanzahl; ++z)
            pipe(ergebnispipe[z]);

        int pid;
        for (int t = 0; t < Threadanzahl; ++t) {
            pid = fork();
            if (pid == 0) {
                MT.seed(time(NULL) + t);
                double erg = 0;
                double ** X = DoubleFeld(N, D);
                int durchlaeufe = 1000000 / Threadanzahl;
                for (int m = 0; m < durchlaeufe; ++m) {
                    double S[1];
                    S[0] = start;
                    Pfadgenerieren(X, 0, S);
                    double sum = 0;
                    for (int n = 0; n < N; ++n)
                        sum += this->semi_f(n, X[n]);
                    erg += sum / (double) durchlaeufe;
                }
                InPipeSchreiben(ergebnispipe[t], erg);
                exit(0);
            }
        }
        double erg = 0;
        double ergebnisse[Threadanzahl];
        for (int f = 0; f < Threadanzahl; ++f) {
            ergebnisse[f] = AusPipeLesen(ergebnispipe[f]);
            erg += ergebnisse[f] / (double) (Threadanzahl);
            if (verbose)printf("Ergebnis %d: %f\n", f, ergebnisse[f]);
        }
        double s = 0;
        for (int f = 0; f < Threadanzahl; ++f)
            s += pow(ergebnisse[f] - erg, 2);
        printf("geschaetzte Varianz: %f\n", s / double(Threadanzahl - 1));
        printf("Testing Ergebnis: %f\n", erg);
        f2 << erg << endl;
    }
    f2.close();*/
}

