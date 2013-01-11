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


AmericanOption* zeiger3;


void* DELEGATE_stuetzerwartung_ausrechnen_THREAD(void* data) {
	zeiger3->stuetzerwartung_ausrechnenThread(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void* DELEGATE_inner_paths_erzeugen_THREAD(void* data) {
	zeiger3->inner_paths_erzeugen_THREAD(((int*)data)[0]);
	pthread_exit(NULL);
	return NULL;
}

void AmericanOption::inner_paths_erzeugen_THREAD(int threadnummer){
	RNG generator;
	double** wdiff=DoubleFeld(N,D);
	int lauf=0;
	for(int u=0;u<durchlaeufe;++u)
		for (int j = 0; j < J; ++j)
			for (int m = 0; m < M; ++m)
				if(m%Threadanzahl==threadnummer)
				{
					lauf++;

					for(int n=0;n<N;++n)
						for(int d=0;d<D;++d)
							if(lauf%2==1)
								wdiff[n][d]=sqrt(dt)*generator.nextGaussian();
							else wdiff[n][d]*=-1.;
					Pfadgenerieren(semi_inner_paths[u][j][m],wdiff, 0,N, stuetzpunkte[j]);
				}
	deleteDoubleFeld(wdiff,N,D);
}



void AmericanOption::stuetzerwartung_ausrechnen(){
	time_t time1 = time(NULL);
	pthread_t threads[Threadanzahl];
	int* nummern=IntFeld(Threadanzahl);
	for (int t = 0; t < Threadanzahl; t++){
		nummern[t]=t;
		pthread_create(&threads[t], NULL, DELEGATE_stuetzerwartung_ausrechnen_THREAD,  &(nummern[t]));
	}
	for (int t = 0; t < Threadanzahl; t++)
		pthread_join(threads[t], NULL);
	deleteIntFeld(nummern,Threadanzahl);
	if(verbose)printf("--------------Time forstuetzerwartung_ausrechnen:%ld seconds\n", time(NULL) - time1);
}

void AmericanOption::semi_inner_paths_erzeugen(){
	if (verbose)printf("innere Pfade erzeugen\n");
	pthread_t threads[Threadanzahl];
	int* nummern=IntFeld(Threadanzahl);
	for (int j = 0; j < Threadanzahl; j++)
	{
		nummern[j]=j;
		pthread_create(&threads[j], NULL, DELEGATE_inner_paths_erzeugen_THREAD,  &(nummern[j]));
	}
	for (int j = 0; j < Threadanzahl; j++)
		pthread_join(threads[j], NULL);
	deleteIntFeld(nummern,Threadanzahl);
}

void AmericanOption::semi() {
	zeiger3 = this;
	mlsm = true;

	Daten(); //Problemdaten laden

	neueExerciseDates(N);
	N=Testing_Dates;  // fuer longstaff schwartz keine hoehere genauigkeit noetig, da explizite formel fuer pfade verwendet wurde
	dt=T/(double(N-1));

	int faktor;
	printf("LSM_k: %d,%d,%d,%d,%d\n",LSM_K0,LSM_K1,LSM_K2,LSM_K3,LSM_K4);

	//    int semi_durchlaeufe=10;   // Wie viele cycles training und testing
	int L=1;
	if (D == 1) {
		Mphi = 56; //16,56
		J = 80; //25,80
		M = 10000; // 5000,10000
		faktor=1;  //1
		L=1;        //Optimierungsversuche 1
		durchlaeufe = 5; //mehrmals pro zeitschritt optimieren 5
		semi_testingpaths = 1e4*10; //Testingpaths, 1e6*10
	}

	if (D == 2) {
		Mphi = 7+3000; //3007   // Basisfunktionen
		J = 200; //200 // Stuetzpunkte
		M = 10000; //10000       // Pfade an jedem stuetzpunkt zum schaetzen
		faktor=2;  //2
		L=10;      //10
		durchlaeufe = 1; //mehrmals pro zeitschritt optimieren 1
		semi_testingpaths = 1e6; // Testingpaths 1e6
	}

	if (D > 2) {
		Mphi = 1+3+D*2+(D-1)+1+2+7000; //+7000
		J = 200; //200
		M = 10000; //10000
		faktor=2;  //2
		L=20; //100
		durchlaeufe = 1;  //1
		semi_testingpaths = 1e5; //Testingpaths 1e6
	}

	printf("Dimensionen: %d\n",D);
	printf("Basisfunktionen: %d \n",Mphi);
	printf("Subsimulation: %d\n",M);
	printf("Stützpunkte: %d\n",J);
	printf("Durchläufe: %d\n",durchlaeufe);

	if (verbose)printf("stuetzpunkte setzen\n");
	stuetzpunkte = DoubleFeld(J, D);
	semi_inner_paths = DoubleFeld(durchlaeufe,J,M,N,D);
	stuetzpunkte_setzen(N/2);
	semi_inner_paths_erzeugen();
	//	stuetzpunkte_ausrichten();

	semi_betas = DoubleFeld(N, Mphi);
	semi_betas_index = IntFeld(N, Mphi);
	semi_betas_index_max = IntFeld(N);
	stuetzerwartung = DoubleFeld(J);

	printf("koeff-testingpfade erzeugen\n");
	RNG generator;
	koeff_testingpaths=DoubleFeld(10000,N,D);
	for(int u=0;u<10000;++u)
		Pfadgenerieren(koeff_testingpaths[u],0,X0,&generator);

	stuetzstelle_active=new bool[J];
	Matrix = DoubleFeld(J, Mphi);
	C = DoubleFeld(Mphi);

	nactual=N-1;
	double training=koeff_testen(NULL);

	//		double* x=DoubleFeld(2);
	//		x[0]=110;x[1]=110;
	//	//	double* grad=semi_Basisfunktionen1DAbl(0,);
	//		double* grad=semi_Basisfunktionen2DAbl(0,900,x);
	//		printf("Test Abl: %f\n",grad[0]);
	//		exit(0);

	for (int n = N - 2; n >= 0; --n) {
		for(int u=0;u<10000;++u)
			Pfadgenerieren(koeff_testingpaths[u],0,X0,&generator); //Achtung, muss nicht sein

		printf("Zeitschritt %d\n", n);
		nactual = n;

		if(D>1){

			//			printf("innere Pfade teilen\n");cout.flush();
			//			if(n!=N-2)
			//				for(int durch=0;durch<durchlaeufe;durch++)
			//				for(int j=0;j<J;++j)
			//					for(int m=0;m<M;++m)
			//						for(int n=0;n<N;++n)
			//							for(int d=0;d<D;++d)
			//								semi_inner_paths[durch][j][m][n][d]/=stuetzpunkte[j][d];

			printf("Stuetzstellen setzen\n");cout.flush();
			stuetzpunkte_setzen(nactual);

			//			printf("innere Pfade multiplizieren\n");cout.flush();
			//			if(n!=N-2)
			//				for(int durch=0;durch<durchlaeufe;durch++)
			//					for(int j=0;j<J;++j)
			//						for(int m=0;m<M;++m)
			//							for(int n=0;n<N;++n)
			//								for(int d=0;d<D;++d)
			//									semi_inner_paths[durch][j][m][n][d]*=stuetzpunkte[j][d];
			//			if(n==N-2)
			semi_inner_paths_erzeugen();
			//		stuetzpunkte_ausrichten();
		}

		double** semi_betas_Feld = DoubleFeld(durchlaeufe, Mphi);
		double** semi_betas_Feld2 = DoubleFeld(durchlaeufe, Mphi);
		for (int i = 0; i < durchlaeufe; ++i) {
			durchlaufactual = i;

			printf("Erwartungen ausrechnen\n");cout.flush();
			stuetzerwartung_ausrechnen();

			if (verbose)stuetzpunkte_ausgeben();

			double min=99999999;
			//			double max=-9999999;
			//double (*pointer)(double*, int);
			//pointer=&(zeiger3->payoff);
			for (int m = 0; m < Mphi; ++m)
				for (int j = 0; j < J; ++j)
					Matrix[j][m] = semi_Basisfunktionen(nactual, m, stuetzpunkte[j]);
			double ** temp_koeff=DoubleFeld(L,Mphi);
			double testergebnisse[L];
			for(lauf=0;lauf<L;++lauf){
				int number_active=0;
				for(int j=0;j<J;++j)
					if(rand()%faktor==0){
						stuetzstelle_active[j]=true;
						number_active++;
					}else
						stuetzstelle_active[j]=false;

				temp_koeff[lauf] =LP_mitGLPK_Loesen(Matrix,true, stuetzerwartung);
				testergebnisse[lauf]=koeff_testen(temp_koeff[lauf]);
				printf("Optimierung (min)%d,\t (%d Stellen aktiv):\t %f\n",lauf,number_active,testergebnisse[lauf]);

				//				double* temp_koeff2 = LP_mitGLPK_Loesen(Matrix,false, stuetzerwartung);
				//				double testergebnis2=koeff_testen(temp_koeff2);
				//				printf("Optimierung (max)%d, (%d Stellen aktiv):\t\t %f\n",lauf,number_active,testergebnis2);
				//				if(testergebnis2>max)
				//				{
				//					deleteDoubleFeld(semi_betas_Feld2[i],Mphi);
				//					semi_betas_Feld2[i] = temp_koeff2;
				//					max=testergebnis2;
				//				}else
				//					deleteDoubleFeld(temp_koeff2,Mphi);
			}

			int* reihe= quicksort((double*)testergebnisse,L);

			for(int m=0;m<Mphi;++m)
				semi_betas_Feld[i][m]=temp_koeff[reihe[0]][m];
			min=testergebnisse[reihe[0]];
			//			if(testergebnis<min)
			//							{
			//								deleteDoubleFeld(semi_betas_Feld[i],Mphi);
			//								semi_betas_Feld[i] = temp_koeff;
			//								min=testergebnis;
			//							}else
			//								deleteDoubleFeld(temp_koeff,Mphi);
			deleteDoubleFeld(temp_koeff,L,Mphi);
			deleteIntFeld(reihe,L);
			printf("Minimum: %f\n",min);

			//			printf("Maximum: %f\n",max);
			training+=min;
		}

		//Durchschnitt als Ergebniss nehmen
		for (int m = 0; m < Mphi; ++m) {
			semi_betas[n][m] = 0;
			for (int i = 0; i < durchlaeufe; ++i)
				semi_betas[n][m] += semi_betas_Feld[i][m] / (double) durchlaeufe;
			//+0.5* semi_betas_Feld2[i][m] / (double) durchlaeufe  ;
		}


		//koeffizienten indizieren
		int indexlauf = 0;
		for (int m = 0; m < Mphi; ++m)
			if (semi_betas[n][m] != 0) {
				semi_betas_index[n][indexlauf] = m;
				indexlauf++;
			}
		semi_betas_index_max[n] = indexlauf;

		printf("Anzahl nichtnegativer Koeff. %d\n", semi_betas_index_max[n]);
		if (verbose)semi_ergebnisse_ausgeben();
		if(n==6 && verbose)exit(0);
		deleteDoubleFeld(semi_betas_Feld,durchlaeufe, Mphi);
		deleteDoubleFeld(semi_betas_Feld2,durchlaeufe, Mphi);
		//exit(0);
	}

	printf("Training %f\n",training);
	semi_testing();
	ErgebnisAnhaengen(training,(char*)"ergebnisse_semi_training.txt");
	deleteDoubleFeld(Matrix,J,Mphi);
	deleteDoubleFeld(C,Mphi);
	deleteDoubleFeld(semi_betas,N, Mphi);
	deleteIntFeld(semi_betas_index ,N, Mphi);
	deleteIntFeld(semi_betas_index_max,N);
	deleteDoubleFeld(stuetzerwartung ,J);
	deleteDoubleFeld(semi_inner_paths,durchlaeufe,J,M,N,D);
	delete[] stuetzstelle_active;
}

double AmericanOption::linearCombinationOfBasis(int zeit, double* x) {
	double sum = 0;
	//    for (int m = 0; m < Mphi; ++m)
	//        sum += semi_betas[zeit][m] * semi_Basisfunktionen(zeit, m, x);
	int m;
	for (int i = 0; i < semi_betas_index_max[zeit]; ++i) {
		m = semi_betas_index[zeit][i];
		sum += semi_betas[zeit][m] * semi_Basisfunktionen(zeit, m, x);
	}
	return sum+0.10*(zeit==7 && verfaelscht);
}

double* AmericanOption::linearCombinationOfBasis_Abl(int zeit, double* x) {
	double* sum =DoubleFeld(D);
	//    for (int m = 0; m < Mphi; ++m)
	//        sum += semi_betas[zeit][m] * semi_Basisfunktionen(zeit, m, x);
	int m;
	for (int i = 0; i < semi_betas_index_max[zeit]; ++i) {
		m = semi_betas_index[zeit][i];
		double* B=semi_BasisfunktionenAbl(zeit, m, x);
		for(int d=0;d<D;++d)
			sum[d] += semi_betas[zeit][m] * B[d];
		deleteDoubleFeld(B,D);
	}
	return sum;
}

double AmericanOption::linearCombination(double* koeff, double* x) {
	double sum = 0;
	for (int m = 0; m < Mphi; ++m)
		sum += koeff[m] * semi_Basisfunktionen(nactual, m, x);
	return sum;
}


double AmericanOption::semi_f(int zeit, double* x) {
	if (zeit == N - 1)
		return payoff(x, zeit);
	return max(0, payoff(x, zeit) - linearCombinationOfBasis(zeit, x));
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


void AmericanOption::stuetzpunkte_setzen(int n) {
	if (D == 1) {
		for (int i = 0; i < J; ++i)
			if (option == MIN_PUT)
				stuetzpunkte[i][0] = Strike * (double) (0.1 + (i + 1) / (double) J);
			else
				stuetzpunkte[i][0] = Strike * (double) (0.9 + (i + 1) / (double) J);
	}

	//	if (D == 332) {
	//		int WJ = (int) (sqrt(J));
	//		for (int i = 0; i < WJ; ++i)
	//			for (int k = 0; k < WJ; ++k) {
	//				stuetzpunkte[i * WJ + k][0] = Strike * (0.01 + 2 * (double) (i) / (double) (WJ));
	//				stuetzpunkte[i * WJ + k][1] = Strike * (0.01 + 2 * (double) (k) / (double) (WJ));
	//			}
	//	}

	//		if (D == 3) {
	//			int WJ = (int) (ceil(pow(J, 1. / 3.)));
	//
	//			for (int i = 0; i < WJ; ++i)
	//				for (int k = 0; k < WJ; ++k)
	//					for (int j = 0; j < WJ; ++j) {
	//						stuetzpunkte[i * WJ * WJ + k * WJ + j][0] = Strike * (0.01 + 2 * (double) (i) / (double) (WJ));
	//						stuetzpunkte[i * WJ * WJ + k * WJ + j][1] = Strike * (0.01 + 2 * (double) (k) / (double) (WJ));
	//						stuetzpunkte[i * WJ * WJ + k * WJ + j][2] = Strike * (0.01 + 2 * (double) (j) / (double) (WJ));
	//					}
	//		}
	//	bool gefaechert=false;

	if (D >1) {
		printf("zufaellige Stuetzstellen\n");cout.flush();
		RNG generator;
		double startwert[D];
		double** X = DoubleFeld(N, D);
		for (int j = 0; j < J; ++j) {
			for(int d=0;d<D;++d)
				startwert[d]=X0[d];//*exp((gefaechert==true)*0.1*generator.nextGaussian());
			Pfadgenerieren(X,0,startwert,&generator);
			for (int d = 0; d < D; ++d)
				stuetzpunkte[j][d] = X[n][d];
		}
		deleteDoubleFeld(X,N,D);
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
		//	for (int j = 0; j < J; ++j)
		//	printf("%.3lf, ", RS[j]);
		printf("\n");
		printf("Funktional Koeff. ausgeben\n");
		for (int j = 0; j < Mphi; ++j)
			printf("%.3lf, ", C[j]);
		printf("\n");
	}
}


double* AmericanOption::LP_mitGLPK_Loesen(double** Matrix, bool kleinergleich, double* RS) {
	time_t time1 = time(NULL);
	//	for (int m = 0; m < Mphi; ++m)
	//		for (int j = 0; j < J; ++j)
	//			if(stuetzstelle_active[j])
	//				Matrix[j][m] = semi_Basisfunktionen(nactual, m, stuetzpunkte[j]);
	//			else
	//				Matrix[j][m]=0;
	//
	double* C=new double[Mphi];

	for (int m = 0; m < Mphi; ++m) {
		C[m] = 0;
		for (int j = 0; j < J; ++j)
			if(stuetzstelle_active[j])
				C[m] += semi_Basisfunktionen(nactual, m, stuetzpunkte[j]);
		if(!kleinergleich)C[m]*=-1;
	}

	if(verbose)printf("Optimization aufstellen time:%ld seconds\n", time(NULL) - time1);

	if (!verbose)glp_term_out(GLP_OFF);
	glp_prob *lp;

	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, J);

	for (int j = 0; j < J; ++j){
		if(kleinergleich)
			glp_set_row_bnds(lp, j + 1, GLP_UP, 0.0, RS[j]);
		else
			glp_set_row_bnds(lp, j + 1, GLP_LO, RS[j],10000000);
	}

	//Zielfunktional uebergeben
	glp_add_cols(lp, Mphi);
	for (int m = 0; m < Mphi; ++m) {
		glp_set_col_bnds(lp, m + 1, GLP_FR, 0.0, 0.0);
		glp_set_obj_coef(lp, m + 1, C[m]);
	}

	//Matrixeintraege uebergeben
	//printf("indexlaenge %d\n",indexlaenge);
	int* ia=new int[1 + J * Mphi];
	int* ja=new int[1 + J * Mphi];
	double* ar=new double[1 + J * Mphi];

	int zaehler = 1;
	for (int i = 0; i < J; ++i)
		for (int j = 0; j < Mphi; j++) {
			ia[zaehler] = i + 1;
			ja[zaehler] = j + 1;
			ar[zaehler] = (stuetzstelle_active[i]==true)*Matrix[i][j];
			zaehler++;
		}

	glp_load_matrix(lp, Mphi*J, ia, ja, ar);
	glp_simplex(lp, NULL);

	//Loesung auslesen
	double* x = DoubleFeld(Mphi);
	for (int m = 0; m < Mphi; ++m)x[m] = glp_get_col_prim(lp, m + 1);

	glp_delete_prob(lp);
	if(verbose)printf("Optimization time:%ld seconds\n", time(NULL) - time1);
	delete[] ja;
	delete[] ia;
	delete[] ar;
	return x;
}

