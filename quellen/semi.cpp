#include "AmericanOption.h"
#include <ctype.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Hilfsmittel.h"
#include <string.h>
#include <glpk.h>

using namespace std;
//using namespace alglib;

AmericanOption* zeiger3;

void* DELEGATE_stuetzerwartung_ausrechnen(void* data) {
	zeiger3->stuetzerwartung_ausrechnenThread(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void* DELEGATE_inner_paths_erzeugen(void* data) {
	zeiger3->inner_paths_erzeugenThread(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void AmericanOption::inner_paths_erzeugenThread(int threadnummer){
	RNG generator;
	for(int u=0;u<durchlaeufe;++u)
		for (int j = 0; j < J; ++j)
			for (int m = 0; m < M; ++m)
				if(m%Threadanzahl==threadnummer)
					Pfadgenerieren(semi_inner_paths[u][j][m], 0, stuetzpunkte[j],&generator);
}

void AmericanOption::semi() {
	zeiger3 = this;
	mlsm = true;
	LSM_setting();

	printf("LSM_k: %d,%d,%d,%d,%d\n",LSM_K0,LSM_K1,LSM_K2,LSM_K3,LSM_K4);


	semi_testingpaths = 1000000; //Wie viele Testingpaths
	//    int semi_durchlaeufe=10;   // Wie viele cycles training und testing


	if (D == 1) {
		Mphi = 1; //56
		J = 80; //80
		M = 10000;
		durchlaeufe = 10; //mehrmals pro zeitschritt optimieren 5
	}

	if (D == 2) {
		Mphi = 1712; //37   // Basisfunktionen
		J = 121; //10*10;//49 // Stuetzpunkte
		M = 500; //5000       // Pfade an jedem stuetzpunkt zum schaetzen
		durchlaeufe = 3; //mehrmals pro zeitschritt optimieren 10
	}

	if (D == 3) {
		Mphi = 1+3+D*2+(D>2?D-1:0)+1+/*8*D*/+1500;
		J = 125; //216   125
		M = 50;   //5000
		durchlaeufe = 5; //mehrmals pro zeitschritt optimieren 5
	}

	if (D > 3) {
		printf("Error 45678\n");
	}

	printf("Dimensionen: %d\n",D);
	printf("Basisfunktionen: %d \n",Mphi);
	printf("Subsimulation: %d\n",M);
	printf("Stützpunkte: %d\n",J);
	printf("Durchläufe: %d\n",durchlaeufe);

	if (verbose)printf("stuetzpunkte setzen\n");
	stuetzpunkte = DoubleFeld(J, D);
	stuetzpunkte_setzen();

	if (verbose)printf("innere Pfade erzeugen\n");

	semi_inner_paths = DoubleFeld(durchlaeufe,J,M,N,D);

	pthread_t threads[Threadanzahl];
	for (int j = 0; j < Threadanzahl; j++)
		pthread_create(&threads[j], NULL, DELEGATE_inner_paths_erzeugen, array_machen(j));
	for (int j = 0; j < Threadanzahl; j++)
		pthread_join(threads[j], NULL);

	semi_betas = DoubleFeld(N, Mphi);
	semi_betas_index = IntFeld(N, Mphi);
	semi_betas_index_max = IntFeld(N);
	stuetzerwartung = DoubleFeld(J);

	for (int n = N - 2; n >= 0; --n) {
		printf("Zeitschritt %d\n", n);
		nactual = n;

		double** semi_betas_Feld = DoubleFeld(durchlaeufe, Mphi);
		for (int i = 0; i < durchlaeufe; ++i) {
			durchlaufactual = i;

			//Stuetzerwartung ausrechnen
			time_t time1 = time(NULL);
			pthread_t threads[Threadanzahl];
			for (int j = 0; j < Threadanzahl; j++)
				pthread_create(&threads[j], NULL, DELEGATE_stuetzerwartung_ausrechnen, array_machen(j));
			for (int j = 0; j < Threadanzahl; j++)
				pthread_join(threads[j], NULL);
			if (verbose)stuetzpunkte_ausgeben();
			time_t time2 = time(NULL);
			if (verbose)printf("Time for Estimation time:%ld seconds\n", time2 - time1);

			//LP aufstellen
			Matrix = DoubleFeld(J, Mphi);
			C = DoubleFeld(Mphi);
			RS = stuetzerwartung;
			semi_betas_Feld[i] = LP_mitGLPK_Loesen(0, J);
			deleteDoubleFeld(Matrix,J,Mphi);
			deleteDoubleFeld(C,Mphi);
		}

		//Durchschnitt als Ergebniss nehmen
		for (int m = 0; m < Mphi; ++m) {
			semi_betas[n][m] = 0;
			for (int i = 0; i < durchlaeufe; ++i)
				semi_betas[n][m] += semi_betas_Feld[i][m] / (double) durchlaeufe;
		}

		int indexlauf = 0;
		for (int m = 0; m < Mphi; ++m)
			if (semi_betas[n][m] != 0) {
				semi_betas_index[n][indexlauf] = m;
				indexlauf++;
			}

		semi_betas_index_max[n] = indexlauf;

		if (verbose)printf("Anzahl nichtnegativer Koeff. %d\n", semi_betas_index_max[n]);
		if (verbose)semi_ergebnisse_ausgeben();

deleteDoubleFeld(semi_betas_Feld,durchlaeufe, Mphi);
	}

	semi_testing();
	deleteDoubleFeld(semi_betas,N, Mphi);
	deleteIntFeld(semi_betas_index ,N, Mphi);
	deleteIntFeld(semi_betas_index_max,N);
	deleteDoubleFeld(stuetzerwartung ,J);
	deleteDoubleFeld(semi_inner_paths,durchlaeufe,J,M,N,D);
}

double AmericanOption::linearCombinationOfBasis(int zeit, double* x, int d) {
	double sum = 0;
	//    for (int m = 0; m < Mphi; ++m)
	//        sum += semi_betas[zeit][m] * semi_Basisfunktionen(zeit, m, x);
	int m;
	for (int i = 0; i < semi_betas_index_max[zeit]; ++i) {
		m = semi_betas_index[zeit][i];
		sum += semi_betas[zeit][m] * semi_Basisfunktionen(zeit, m, x);
	}
	return sum;
}

double Min(double x, double y) {
	return x < y ? x : y;
}

double Min(double x, double y, double z) {
	return min(x, min(y, z));
}

double AmericanOption::semi_f(int zeit, double* x) {
	if (zeit == N - 1)
		return payoff(x, zeit);
	//----------------------------
	//	if(nactual==N-1)
	int dmax = argMax(x, D);
	dmax = 0;
	return max(0, payoff(x, zeit) - linearCombinationOfBasis(zeit, x, dmax));

	//return max(0, payoff(x,zeit) - linearCombinationOfBasis(zeit,x) );
	/*
    //	else
            -------------
            //			return Min(
            //					max(payoff(x,zeit) - european_MaxCall_ND(x,(double)zeit*dt,(double)zeit*dt+dt),0),
            //					max(payoff(x,zeit) - european_MaxCall_ND(x,(double)zeit*dt,T)                 ,0),
            //					max(payoff(x,zeit) - linearCombinationOfBasis(zeit,x)                         ,0)
            //			);
            //----------------------------
            //	return max(0, payoff(x,zeit) -european_MaxCall_ND(x,zeit*dt,zeit*dt+dt));
            //printf("unterschied: %f\n",linearCombinationOfBasis(zeit,x) - LSM_C_estimated(x,zeit));
            //	return max(0, payoff(x,zeit) - LSM_C_estimated(x,zeit)); //TODO
            //		return max(0, payoff(x,zeit) - E);
            //return max(0, linearCombinationOfBasis(zeit,semi_betas[zeit],x));
	 */
}

void AmericanOption::stuetzerwartung_ausrechnenThread(int k) {
	for (int j = 0; j < J; ++j)
		if (j % Threadanzahl == k) { // jeder Thread muss nur einige bearbeiten
			stuetzerwartung[j] = 0;
			for (int m = 0; m < M; ++m)
				for (int nnn = nactual + 1; nnn < N; ++nnn)
					stuetzerwartung[j] += semi_f(nnn, semi_inner_paths[durchlaufactual][j][m][nnn - nactual]) / (double) (M);
		}
}

void AmericanOption::stuetzpunkte_setzen() {

	if (D == 1) {
		for (int i = 0; i < J; ++i)
			if (option == MIN_PUT)
				stuetzpunkte[i][0] = Strike * (double) (0.1 + (i + 1) / (double) J);
			else
				stuetzpunkte[i][0] = Strike * (double) (0.9 + (i + 1) / (double) J);
	}

	if (D == 2) {
		int WJ = (int) (sqrt(J));
		//int WJ=7;
		for (int i = 0; i < WJ; ++i)
			for (int k = 0; k < WJ; ++k) {
				stuetzpunkte[i * WJ + k][0] = Strike * (0.01 + 2 * (double) (i) / (double) (WJ));
				stuetzpunkte[i * WJ + k][1] = Strike * (0.01 + 2 * (double) (k) / (double) (WJ));
			}
		//stuetzpunkte[30][0]=110;
		//stuetzpunkte[30][1]=60;
	}

	if (D == 3) {
		int WJ = (int) (ceil(pow(J, 1. / 3.)));
		//            printf("sfd %d\n",WJ);exit(0);
		for (int i = 0; i < WJ; ++i)
			for (int k = 0; k < WJ; ++k)
				for (int j = 0; j < WJ; ++j) {
					stuetzpunkte[i * WJ * WJ + k * WJ + j][0] = Strike * (0.01 + 2 * (double) (i) / (double) (WJ));
					stuetzpunkte[i * WJ * WJ + k * WJ + j][1] = Strike * (0.01 + 2 * (double) (k) / (double) (WJ));
					stuetzpunkte[i * WJ * WJ + k * WJ + j][2] = Strike * (0.01 + 2 * (double) (j) / (double) (WJ));
				}
	}

	if (D > 3) {
		printf("Error 4678");
		exit(0);
		//        for (int j = 0; j < J; ++j) {
		//                double** X = DoubleFeld(N, D);
		//            Pfadgenerieren(X);
		//            for (int d = 0; d < D; ++d)
		//                stuetzpunkte[j][d] = X[N / 2][d];
	}
}

void AmericanOption::lp_ausgeben() {
	if (verbose) {
		printf("Matrix ausgeben\n");
		for (int j = 0; j < J; ++j) {
			for (int m = 0; m < Mphi; ++m)
				printf("%.3lf, ", Matrix[j][m]);
			printf("\n");
		}
		printf("RS ausgeben\n");
		for (int j = 0; j < J; ++j)
			printf("%.3lf, ", RS[j]);
		printf("\n");
		printf("Funktional Koeff. ausgeben\n");
		for (int j = 0; j < Mphi; ++j)
			printf("%.3lf, ", C[j]);
		printf("\n");
	}
}

double * semi_ergebnisse;

void* DELEGATE_semi_test(void* data) {
	zeiger3->semi_testThread(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void AmericanOption::semi_testThread(int threadnummer) {
	semi_ergebnisse[threadnummer] = 0;

	double** x = DoubleFeld(N, D);

	int seed = time(NULL) + threadnummer + getpid();
	srand(seed);

	int durchlaeufe = (int)(double)(semi_testingpaths) / (double)(Threadanzahl);
	RNG generator;
	for (int m = 0; m < durchlaeufe; ++m) {
		Pfadgenerieren(x, 0,X0,&generator);
		double sum = 0;
		for (int n = 0; n < N; ++n)
			sum += semi_f(n, x[n]);
		semi_ergebnisse[threadnummer]+=sum/(double)(durchlaeufe);

	}
}

void AmericanOption::semi_testing() {
	printf("Testing\n");

	semi_ergebnisse=DoubleFeld(Threadanzahl);

	pthread_t threads[Threadanzahl];
	for (int j = 0; j < Threadanzahl; j++)
		pthread_create(&threads[j], NULL, DELEGATE_semi_test, array_machen(j));
	for (int j = 0; j < Threadanzahl; j++)
		pthread_join(threads[j], NULL);

	double erg=0;
	for(int threadnummer=0;threadnummer<Threadanzahl;++threadnummer)
		erg+=semi_ergebnisse[threadnummer]/(double)Threadanzahl;

	printf("%f\n", erg);
	ErgebnisAnhaengen(erg,(char*)"ergebnisse_semi.txt");

	//---------------
	/*
    for (int testlauf = 0; testlauf < 10; ++testlauf) {
        int ergebnispipe[Threadanzahl][2];
        for (int z = 0; z < Threadanzahl; ++z)
            pipe(ergebnispipe[z]);

        int pid;
        for (int t = 0; t < Threadanzahl; ++t) {
            pid = fork();
            if (pid == 0) {
                double** x = DoubleFeld(N, D);
                double** wdiff = DoubleFeld(N, D);
                double** sprue = DoubleFeld(N, D);

                int seed = time(NULL) + t + MT() + getpid();
                if (testlauf > 0)seed += (int) (testergs[testlauf - 1]*10000);
                MT.seed(seed);
                double erg = 0;
                int durchlaeufe = 100000 / Threadanzahl; //Achtung
                for (int m = 0; m < durchlaeufe; ++m) {
                    for (int n = 0; n < N; ++n)
                        for (int j = 0; j < D; ++j) {
                            if (m % 2 == 0)wdiff[n][j] = sqrt(dt) * nextGaussian();
                            else wdiff[n][j] = -wdiff[n][j]; //antithetics
                            sprue[n][j] = 0;
                        }
                    Pfadgenerieren(x, wdiff, sprue);
                    double sum = 0;
                    for (int n = 0; n < N; ++n)
                        if (Mphi > 1) {
                            sum += semi_f(n, x[n]);
                        } else { //if(Mphi==1)
                            if (n == N - 1)
                                sum += max(payoff(x[n], n), 0);
                            else {
                                double d[5];
                                d[0] = europeanValue(x[n], (double) n*dt, (double) n * dt + dt);
                                d[1] = europeanValue(x[n], (double) n*dt, T);
//                                d[2] = europeanValue(x[n], (double) n*dt, T * 0.8 + (double) n * dt * 0.2);
//                                d[3] = europeanValue(x[n], (double) n*dt, T * 0.6 + (double) n * dt * 0.3);
//                                d[4] = europeanValue(x[n], (double) n*dt, T * 0.3 + (double) n * dt * 0.7);
//                                sum += max(payoff(x[n], n) - max(max(max(max(d[0], d[1]), d[2]), d[3]), d[4]), 0);
                            sum += max(payoff(x[n], n) - max(d[0], d[1]),0);
                            }
                        }
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
        }
        testergs[testlauf] = erg;
        double zwischen = 0;
        for (int i = 0; i <= testlauf; ++i)
            zwischen += testergs[i] / double(testlauf + 1);
        printf("Zwischenergebnis %d: %f (bei %d Pfaden)\r", testlauf, zwischen, 100000 * (testlauf + 1));

    }
    double erg = 0;
    for (int i = 0; i < 10; ++i)
        erg += testergs[i]*0.1;*/
	//
	////    printf("\n\nGesamtergebnis: %f\n", erg);
	//    ErgebnisAnhaengen(erg);
}

//void AmericanOption::stuetzerwartungen_ausrechnen(){
////	pthread_t threads[J];
////
////	for(int  i = 0; i < J; i++ ) {
////		int* z=(int*)malloc(sizeof(int));z[0]=i;
////		pthread_create( &threads[i], NULL, DELEGATE_stuetzerwartung_ausrechnen, (void*)z );
////		//                       Argument für thread_function ---^
////	}
////
////	// Hier wird gewartet bis alle Threads beendet sind.
////	for(int  i = 0; i < J; i++ )
////		pthread_join(threads[i], NULL);
//
//
//	return;
//	//***************************************************************************
//	//	Threadanzahl=J;
//	//	int ergPipe[Threadanzahl][2];
//	//
//	//	for(int z=0;z<Threadanzahl;++z)
//	//		pipe(ergPipe[z]);
//	//
//	//	int pid;
//	//	for(int t=0;t<Threadanzahl;++t)
//	//	{
//	//		pid = fork();  // Prozess duplizieren
//	//		if (pid == 0)
//	//		{
//	//			double stue=0;
//	//			srand(getpid()+t+time(NULL));
//	//			MT.seed(getpid()+t+time(NULL));
//	//
//	//			double** x    =DoubleFeld(N,D);
//	//			double** wdiff=DoubleFeld(N,D);
//	//			double** sprue=DoubleFeld(N,D);
//	//			stue=0;
//	//			int durchlaeufe=20000;
//	//			if(payoff(stuetzpunkte[t],nactual)==0)durchlaeufe/=10;
//	//			for(int m=0;m<durchlaeufe;++m)
//	//			{
//	//				for(int n=0;n<N;++n)
//	//					for(int j=0;j<D;++j){
//	//						if(m%2==0)wdiff[n][j]=sqrt(dt)*nextGaussian();
//	//						else wdiff[n][j]=-wdiff[n][j]; //antithetics
//	//						sprue[n][j]=0;
//	//					}
//	//				Pfadgenerieren(x,wdiff,sprue,nactual,stuetzpunkte[t]);
//	//				for(int nnn=nactual+1;nnn<N;++nnn)
//	//					stue+=semi_f(nnn,x[nnn])/(double)(durchlaeufe);
//	//			}
//	//			InPipeSchreiben(ergPipe[t],stue);
//	//			exit(0);
//	//		}
//	//	}
//	//	for(int t=0;t<Threadanzahl;++t)  //Ergebnisse aus pipe auslesen
//	//		stuetzerwartung[t]=AusPipeLesen(ergPipe[t]);
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

double* AmericanOption::LP_mitGLPK_Loesen(int* index, int indexlaenge) {
	for (int m = 0; m < Mphi; ++m)
		for (int j = 0; j < J; ++j)
			Matrix[j][m] = semi_Basisfunktionen(nactual, m, stuetzpunkte[j]);

	for (int m = 0; m < Mphi; ++m) {
		C[m] = 0;
		for (int j = 0; j < J; ++j)
			C[m] += semi_Basisfunktionen(nactual, m, stuetzpunkte[j]);
	}



	if(index==NULL){
		indexlaenge=Mphi;
		index=new int[Mphi];
		for(int m=0;m<Mphi;++m)
			index[m]=m;
	}

	if (!verbose)glp_term_out(GLP_OFF);
	glp_prob *lp;

	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, J);

	for (int j = 0; j < J; ++j)
		glp_set_row_bnds(lp, j + 1, GLP_UP, 0.0, RS[j]);

	//Zielfunktional uebergeben
	glp_add_cols(lp, indexlaenge);
	for (int m = 0; m < indexlaenge; ++m) {
		glp_set_col_bnds(lp, m + 1, GLP_FR, 0.0, 0.0);
		glp_set_obj_coef(lp, m + 1, C[index[m]]);
	}

	//Matrixeintraege uebergeben
	int ia[1 + J * indexlaenge], ja[1 + J * indexlaenge];
	double ar[1 + J * indexlaenge];
	int zaehler = 1;
	for (int i = 0; i < J; ++i)
		for (int j = 0; j < indexlaenge; j++) {
			ia[zaehler] = i + 1;
			ja[zaehler] = j + 1;
			ar[zaehler] = Matrix[i][index[j]];
			zaehler++;
		}

	glp_load_matrix(lp, indexlaenge*J, ia, ja, ar);
	glp_simplex(lp, NULL);

	//Loesung auslesen
	double* x = DoubleFeld(indexlaenge);
	for (int m = 0; m < indexlaenge; ++m)x[m] = glp_get_col_prim(lp, m + 1);

	glp_delete_prob(lp);
	return x;
}

void AmericanOption::semi_ergebnisse_ausgeben(){
	printf("semi_betas ausgeben \n");
	for (int j = 0; j < Mphi; ++j)
		printf("%.4lf, ", semi_betas[nactual][j]);
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
			f << linearCombinationOfBasis(nactual, point, 0) << endl;
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
		array_Q[j]= max(payoff(stuetzpunkte[j], nactual) - linearCombinationOfBasis(nactual, stuetzpunkte[j], 0), 0) + stuetzerwartung[j];
		array_payoff[j] =payoff(stuetzpunkte[j], nactual);
		array_diff[j] =linearCombinationOfBasis(nactual, stuetzpunkte[j], 0) - stuetzerwartung[j] ;
		array_stue[j] =stuetzerwartung[j];
		array_linComb[j] = linearCombinationOfBasis(nactual, stuetzpunkte[j], 0);
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

//double AmericanOption::phiPOS(double* x, int j, int time)
//{
//	//return phi(x,j,time);
//	//return max(0,semi_Basisfunktionen(x,j,time));
//}

//		double x0[D];
//		for(int d=0;d<D;++d)
//			x0[d]=K;
//		for (int m = 0; m < M; ++m)
//			Pfadgenerieren(X[m],n,x0);

//		//Stoch Opt
//		double opt=10000000;
//		for(int k=0;k<Mphi;k++)
//			semi_betas[n][k]=0;
//
//		for(int lauf=0;lauf<1000;++lauf)
//		{
//			double gauss[Mphi];
//			for(int i=0;i<Mphi;++i){
//				gauss[i]=0.01*nextGaussian();
//				if(gauss[i]<semi_betas[n][i])gauss[i]=0;
//			}
//
//			for(int k=0;k<Mphi;k++)
//				semi_betas[n][k]+=gauss[k];
//
//			if(verbose){
//				printf("semi_betas=");
//				for(int i=0;i<Mphi;++i)
//					printf("%.2lf, ",semi_betas[n][i]);
//				printf("\n");
//			}
//
//			bool ueber=true;
//			for(int j=0;j<J;++j)
//				if( linComb(n,semi_betas[n],stuetzpunkte[j])+stuetzerwartung[j]<payoff(stuetzpunkte[j],n) )ueber=false;
//			if(verbose)
//				printf("ueber: %d, opt=%f\n",ueber,opt);
//			double wert=linComb(n,semi_betas[n],x0);
//			if(wert<opt && ueber)
//				opt=wert;
//			else
//				for(int k=0;k<Mphi;k++)
//					semi_betas[n][k]-=gauss[k];
//		}


//		for(int j=0;j<J;++j){
//			double ** x=DoubleFeld(N,D);
//			double** wdiff=DoubleFeld(N,D);
//			double** sprue=DoubleFeld(N,D);
//			stuetzerwartung[j]=0;
//			int durchlaeufe=30000;
//			for(int m=0;m<durchlaeufe;++m)
//			{
//				for(int n=0;n<N;++n)
//					for(int j=0;j<D;++j){
//						if(m%2==0)wdiff[n][j]=sqrt(dt)*nextGaussian();
//						else wdiff[n][j]=-wdiff[n][j];
//						sprue[n][j]=0;
//					}
//				Pfadgenerieren(x,wdiff,sprue,nactual,stuetzpunkte[j]);
//				for(int nnn=nactual+1;nnn<N;++nnn)
//					stuetzerwartung[j]+=f(nnn,x[nnn])/(double)(durchlaeufe);
//
//			}
//		}


//bool fertig=false;
////			while(!fertig)
//			{
////				for(int d=0;d<D;++d)
////				stuetzpunkte[t][]=X0[d]*exp( ( (r-delta)-0.5*sigma[0]*sigma[0])*dt*nactual+sqrt(dt*nactual)*nextGaussian());
//			for(int m=0;m<durchlaeufe;++m)
//			{
//				if(m%1000==0){
//					if(stue*durchlaeufe/m<0.001){
//						break;
//						fertig=true;
//					}
//					if(stue*durchlaeufe/m>1.1*payoff(stuetzpunkte[t],nactual) )break;
//				}
//
//


//	//nur zum ueberpruefen
//	string line;
//	ifstream myfile ("ERG.dat");
//	int k=0;
//	if (myfile.is_open())
//	{
//		while ((! myfile.eof())&& k<Mphi+1  )
//		{
//			//printf("gelesen\n");
//			getline (myfile,line);
//			char * buffer = new char[line.length()];
//			strcpy(buffer,line.c_str());
//			//			printf("Empfangene Daten: %s\n",buffer);
//			if(k>0)semi_betas[N-1][k-1]=(double)(atof(buffer));
//			k++;
//		}
//		myfile.close();
//	}
//	else printf("Kann Datei nicht oeffnen\n");


//double AmericanOption::objfsO(double* koeff){
//	double erg=0;
//	for(int j=0;j<J;++j)
//		erg+=-linearCombinationOfBasis(nactual,koeff,stuetzpunkte[j]);
//		//	int j=3;
//		//				erg+=max(payoff(stuetzpunkte[j],nactual)-linearCombinationOfBasis(nactual,koeff,stuetzpunkte[j]),0);
////		erg+=payoff(stuetzpunkte[j],nactual)-linearCombinationOfBasis(nactual,koeff,stuetzpunkte[j]);
//	return erg;
//}



//double* AmericanOption::LP_mitR_Loesen(){
//	//exit(0);
//	double* e=(double*)malloc(sizeof(double)*(Mphi));
//
//	//checken, ob NullMatrix vorliegt
//	bool nichtnull=false;
//	for(int j=0;j<J;++j)
//		for(int m=0;m<Mphi;++m)
//			if(Matrix[j][m]!=0)nichtnull=true;
//	if(!nichtnull)
//	{
//		for(int m=0;m<Mphi;++m)
//			e[m]=0;
//		return e;
//	}
//
//	fstream f;
//	f.open("Matrix.dat", ios::out);
//	for(int j=0;j<J;++j)	{
//		for(int m=0;m<Mphi;++m)
//			f<<Matrix[j][m]<<" ";
//		f<<endl;
//	}
//	f.close();
//
//	f.open("RS.dat", ios::out);
//	for(int j=0;j<J;++j)
//		f<<RS[j]<<" ";
//	f<<endl;
//	f.close();
//
//	f.open("C.dat", ios::out);
//	for(int j=0;j<Mphi;++j)
//		f<<C[j]<<" ";
//	f<<endl;
//	f.close();
//
//	int pid;
//	int status;
//	pid = fork();
//	if (pid == 0){
//		char* path =(char*) malloc(sizeof(char)*127);
//		path = getenv("PWD");
//		if(path!=NULL) {
//			strcpy(path+(strlen(path)),"/lp.sh");
//		}
//		else {
//			printf("Couldn't find script\n");
//			exit(1);
//		}
//		//printf("Executing script : %s\n",path);
//		if(system(path)){
//			perror("error\n");
//		}
//		exit(0);
//	}
//	if (pid > 0) {
//		//  printf("Ich bin der Elternprozess, das Kind ist %d.\n",pid);
//		pid = wait(&status);
//		//		printf("Ende des Prozesses %d: ", pid);
//		if (WIFEXITED(status)) {
//			//  printf("Der Prozess wurde mit exit(%d) beendet.\n"WEXITSTATUS(status));
//		}
//		if (WIFSIGNALED(status)) {
//			printf("Der Prozess wurde mit kill -%d beendet.\n",
//					WTERMSIG(status));
//		}
//	}
//	usleep(50000);
//	bool fertig=false;
//	while(!fertig){
//		usleep(10000);
//		string line;
//		ifstream myfile ("ERG.dat");
//		int k=0;
//		if (myfile.is_open())
//		{
//			while ((! myfile.eof())&& k<Mphi+1  )
//			{
//				//printf("gelesen\n");
//				getline (myfile,line);
//				char * buffer = new char[line.length()];
//				strcpy(buffer,line.c_str());
//				//				printf("Empfangene Daten: %s\n",buffer);
//				if(k>0)e[k-1]=(double)(atof(buffer));
//				k++;
//			}
//			myfile.close();
//		}
//		else printf("Kann Datei nicht oeffnen\n");
//
//		for(int j=0;j<Mphi;++j)
//			if(e[j]!=semi_betas[nactual+1][j])fertig=true;
//		//if(!fertig)printf("auf Datei warten\n");
//	}
//
//	//if(verbose)for(int k=0;k<Mphi;++k)
//	//	printf("e[%d]=%f, ",k,e[k] );
//
//	return e;
//}

//		//ergebnis testen
//		if(verbose){
//			double fehlersumme=0;
//			for(int j=0;j<J;++j)
//				fehlersumme+=pow(linearCombinationOfBasis(n,stuetzpunkte[j])-stuetzerwartung[j],2);
//			printf("Fehlersumme %f\n", fehlersumme);
//		}

//		double* test=LP_mitR_Loesen();
//		for(int i=0;i<Mphi;++i)
//			printf("unterschied %f\n",semi_betas[n][i]-test[i]);



//double* AmericanOption::LP_mitALGLIB_Loesen(){
//	/*
//	real_2d_array mat;
//	mat.setlength(J,Mphi+1);
//
//	for(int j=0;j<J;++j)    // b setzen
//		//			Matrix[j][Mphi]=max(payoff(stuetzpunkte[j],n)-stuetzerwartung[j],0);
//		mat[j][Mphi]=-RS[j];
//
//	for(int m=0;m<Mphi;++m)		// Matrix setzen
//		for(int j=0;j<J;++j)
//			mat[j][m]=-Matrix[j][m];
//
//	//		if(verbose)
//	//		{
//	//			printf("Matrix ausgeben\n");
//	//			for(int j=0;j<J;++j)
//	//			{
//	//				for(int m=0;m<Mphi;++m)
//	//					printf("%.2lf, ",Matrix[j][m]);
//	//				printf("\n");
//	//			}
//	//		}
//
//	real_1d_array x;x.setlength(Mphi);
//	for(int i=0;i<x.length();++i)
//		x[i]=0;
//	//	x[1]=0;
//
//	//	double* koeff=(double*)malloc(sizeof(double)*x.length());
//	//	for(int i=0;i<x.length();++i)
//	//		koeff[i]=x[i];
//
//	integer_1d_array ct;ct.setlength(J);
//	for(int i=0;i<ct.length();++i)
//		ct[i]=1;
//	minbleicstate state; minbleicreport rep;
//	double epsg = 0.000001;	double epsf = 0; double epsx = 0;
//	double diffstep = 1.0e-6; double epso = 0.00001; double epsi = 0.00001;
//	minbleiccreatef(x,diffstep, state);minbleicsetlc(state, mat, ct);
//	minbleicsetinnercond(state, epsg, epsf, epsx);minbleicsetoutercond(state, epso, epsi);
//	alglib::minbleicoptimize(state, function1_func);minbleicresults(state, x, rep);
//	printf("result (4 = korrekt beendet): %d\n", int(rep.terminationtype)); // EXPECTED: 4
//
//	double* erg=DoubleFeld((Mphi));
//	for(int k=0;k<Mphi;k++)
//		erg[k]=x[k];
//	return erg;*/
//	return NULL;
//}



//		double** pts = faurepts(10, J+10, D, 5);
//		for(int n=0;n<100;++n){
//			for(int d=0;d<D;++d){
//				if(pts[n][d]>1.)printf("ERROR 112\n");
//				 printf("%f ",pts[n][d]);
//			}
//			printf("\n");
//		}

//        if(j<30){
//            double summe=0;
//        double p=0.2*((double)j+1.)/30.;
//        for(int d=0;d<D;++d)
//           summe+=exp(p*EuropeanCall1D_discounted(zeit*dt,zeit*dt+dt,x[d],Strike));
//        return 1./p*log(summe);
//        }
//        j-=30;

//        double xx[2];
//        xx[0] = x[0]*2;
//        xx[1] = x[1]*2;
//        for (int i = 0; i < 30; ++i) {
//            if (j == i)return payoff(xx, zeit);
//            xx[0] *= 0.95;
//            xx[1] *= 0.95;
//        }
//        j -= 30;
//
//        for (int i = 0; i < 50; ++i) {
//            double p = (double) (i + 1) / 50. * 0.1;
//            if (j == i)return 1. / p * log(exp(p * x[0]) + exp(p * x[1]));
//        }
//        j -= 10;
//
//        if (j % 2 == 0)return cos((j / 2)*(10. * fabs(x[0] / X0[0] + x[1] / X0[1])));
//        if (j % 2 == 1)return sin((j - 1) / 2 * (10. * fabs(x[0] / X0[0] + x[1] / X0[1])));
//        j -= 10;
//
//        if (j % 2 == 0)return cos((j / 2)*(10. * fabs(x[0] / X0[0] - x[1] / X0[1])));
//        if (j % 2 == 1)return sin((j - 1) / 2 * (10. * fabs(x[0] / X0[0] - x[1] / X0[1])));
//        j -= 10;


//if(argMax(stuetzpunkte[j],D)==d && payoff(stuetzpunkte[j],nactual)!=0)
//if(argMax(stuetzpunkte[j],D)==d)//if(payoff(stuetzpunkte[j],nactual)!=0)
//C[m]+=(1./0.1+payoff(stuetzpunkte[j],nactual))*semi_Basisfunktionen(n,m,stuetzpunkte[j]);
//if(argMax(stuetzpunkte[j],D)==d && payoff(stuetzpunkte[j],nactual)!=0)
//if(argMax(stuetzpunkte[m],D)==d)
//double gewicht=exp(-100.*fabs(payoff(stuetzpunkte[j],nactual)-stuetzerwartung[j]));
//if(D==1)gewicht=1.;

// int indexlaenge = 1000;
//   time_t time3 = time(NULL);
//semi_betas_Feld[i] = LP_mitGLPK_Loesen(NULL,0);
//   time_t time4 = time(NULL);

//            int* index = new int[indexlaenge];
//            for (int m = 0; m < indexlaenge; ++m)
//                index[m] = LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4 + m;
//
//            if (verbose) {
//                printf("Fuer die erste optimierung:\n");
//                for (int m = 0; m < indexlaenge; ++m)
//                    printf("%d, ", index[m]);
//                printf("\n");
//            }
//

//
//            time_t time3 = time(NULL);
//            semi_betas_Feld[i] = LP_mitGLPK_Loesen(index, indexlaenge);
//            time_t time4 = time(NULL);
//            if (verbose)printf("Time for Optimization:%ld seconds\n", time4 - time3);
//
//            int indexlauf = 0;
//            index = new int[Mphi];
//            for (int m = 0; m < LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4; ++m) {
//                index[indexlauf] = m;
//                indexlauf++;
//            }
//            for (int m = 0; m < 1000; ++m)
//                if (semi_betas_Feld[i][m] != 0) {
//                    index[indexlauf] = m + LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4;
//                    indexlauf++;
//                }
//            indexlaenge = indexlauf;
//            if (verbose) {
//                printf("Fuer die zweite optimierung:\n");
//                for (int m = 0; m < indexlaenge; ++m)
//                    printf("%d, ", index[m]);
//                printf("\n");
//            }
