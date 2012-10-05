
#include "AmericanOption.h"

using namespace std; 

double AmericanOption::semi_Basisfunktionen(int zeit, int j, double* x) {
	double y[D];
	for (int d = 0; d < D; ++d) {
		y[d] = x[d];
		y[d] -= 100;
		y[d] /= 100;
	}

	if (D == 1) {
		if (j < 6) {
			if (j == 0)
				return 1;
			if (j == 1)
				return y[0];
			if (j == 2)
				return y[0] * y[0];
			if (j == 3)
				return y[0] * y[0] * y[0];
			if (j == 4)
				return y[0] * y[0] * y[0] * y[0]; //Polynom 5. Grades bringt nur Probleme
				if (j == 5)
					return payoff(x, zeit);
		}
		j -= 6;

		if (option == MIN_PUT)
			return max(Strike * 4. * (double) (j) / (double) (Mphi - 6) - x[0], 0);
		else
			return max(x[0] - Strike * (double) (j) / (double) (Mphi - 6), 0);
	}

	if (D == 2) {
		//           if (j == 0)
		//               return 1;
		//           if (j == 1)
		//               return europeanValue(x,zeit*dt,zeit*dt+dt);
		//           if (j == 2)
		//               return europeanValue(x,zeit*dt,T);
		//             if (j == 3)
		//                return payoff(x, zeit);
		//           j-=4;

		if (j == 0)
			return 1;
		if (j == 1)
			return y[0];
		if (j == 2)
			return y[1];
		if (j == 3)
			return y[0] * y[1];
		if (j == 4)
			return y[0] * y[0];
		if (j == 5)
			return y[1] * y[1];
		if (j == 6)
			return payoff(x, zeit);
		j -= 7;

		if (j == 0)
			return sqrt(y[0]);
		if (j == 1)
			return sqrt(y[1]);
		if (j == 2)
			return sqrt(y[0] * y[1]);
		if (j == 3)
			return sqrt(y[0] * y[0]);
		if (j == 4)
			return sqrt(y[1] * y[1]);
		j -= 5;

		if (j == 0)
			return exp(-y[0]);
		if (j == 1)
			return exp(-y[1]);
		if (j == 2)
			return exp(-y[0] + y[1]);
		if (j == 3)
			return exp(-y[0] + y[0]);
		if (j == 4)
			return exp(-y[1] + y[1]);
		j -= 5;


		//            if (j == 0)
			//                return pow(y[0]-y[1],1);
		//            if (j == 1)
			//                return pow(y[0]-y[1],2);
		//            if (j == 2)
			//                return pow(y[0]-y[1],3);
		//            if (j == 3)
			//                return pow(y[0]+y[1],1);
		//            if (j == 4)
		//                return pow(y[0]+y[1],2);
		//            if (j == 5)
		//                return pow(y[0]+y[1],3);
		//        j -= 6;

		//  double xx[2];
		//  xx[0] = x[0]*2;
		//  xx[1] = x[1]*2;
		//        for (int i = 0; i < 300; ++i)// {
		//            if (j == i)return payoff(xx, zeit);
		//   xx[0] *= 0.99;
		//    xx[1] *= 0.99;
		//  }
		if (j < 700)
			return max(0.5 * (x[0] + x[1]) - 4 * j / 700. * X0[0], 0);
		j -= 700;

		if (j < 1000) {
			double xx[2];
			//            xx[0]=x[0]*2.*pow(0.99,j);
			//            xx[1]=x[1]*2.*pow(0.99,j);
			xx[0] = j / 1000. * 2. * x[0];
			xx[1] = j / 1000. * 2. * x[1];
			return payoff(xx, zeit);
		}
		j -= 1000;


		//       for (int i = 0; i < 100; ++i) {
		//            double p = (double) (i + 1) / 100. * 0.1;
		//            if (j == i)return 1. / p * log(exp(p * x[0]) + exp(p * x[1]));
		//        }
		//        j -= 100;
	}

	if (D > 2) {
		//if (j < LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4)return LSM_phi(x, j, zeit);
	//	j -= LSM_K0 + LSM_K1 + LSM_K2 + LSM_K3 + LSM_K4;

		//        if(j==0){
		//            return europeanValue(x,zeit*dt,T);
		//        }
		//
		//        if(j==1){
		//            return europeanValue(x,zeit*dt,zeit*dt+dt);
		//        }

		if (j < 400) {
			double xx[D];
			for (int d = 0; d < D; ++d)
				xx[d] = (double) j / 400. * 4. * x[d];
			return payoff(xx, zeit);
		}
		j -= 400;

		if (j < 400) {
			double summe = 0;
			for (int d = 0; d < D; ++d)
				summe += x[d];
			return max(summe - (double) (j) / 400. * 8. * (double) D * X0[0], 0);
		}
		j -= 400;

		if (j < 200) {
			double summe = 0;
			for (int d = 0; d < D; ++d)
				if(d!=j%3)
					summe += x[d];
			return max(summe - (double) (j) / 200. * 6. * (double) D * X0[0], 0);
		}
		j -= 200;

//		if (j < 200) {
//			//double summe = 0;
//			// for (int d = 0; d < D; ++d)
//				//    if(d!=j%3)
//					//   summe += x[d];
//			//return exp (x[0]/100.*(double) (j) / 200.) * exp (x[1]/100.*(double) (j) / 200.);
//			return   exp((x[0]/100.-0.5)*(double) (j) / 200.)
//					*exp((x[1]/100.-0.5)*(double) (j) / 200.)
//					*exp((x[2]/100.-0.5)*(double) (j) / 200.);
//
//		}
//		j -= 200;


//		if (j < 200) {
//			int reihe[D];
//			for(int jj=0;jj<D;++jj)
//				reihe[jj]=jj;
//			BubbleSort(x,reihe,D);
//			return pow(x[reihe[0]]-x[reihe[1]],(double)j/200.*4);
//		}
//		j -= 200;
	}

	printf("Error598");
	return -1;
}


//        if (j == 0)
//            return exp(-y[0]);
//        if (j == 1)
//            return exp(-y[1]);
//        if (j == 2)
//            return exp(-y[2]);
//        if (j == 3)
//            return exp(-y[0] - y[1] - y[2]);
//        if (j == 4)
//            return exp(-y[0] * y[1]);
//        if (j == 5)
//            return exp(-y[1] * y[2]);
//        if (j == 6)
//            return exp(-y[0] * y[1]);
//        j -= 7;
