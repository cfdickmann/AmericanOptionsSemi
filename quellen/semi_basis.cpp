
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
		int reihe[D];
		for(int jj=0;jj<D;++jj)
			reihe[jj]=jj;
		BubbleSort(x,reihe,D);

		if(j<1)return 1;
		j-=1;

		if(j<3)
			return pow(x[reihe[0]],j+3);
		j-=3;

		if(j<D*2 )
		{
			int a=j%D;
			int b=(j-a)/D;
			//if(verbose)printf("%d:asset auf platz %d hoch %d\n",j,a,b);
			return pow(x[reihe[a]],b+1);
		}
		j-=D*2;

		if(j<(D>2?D-1:0))
			return x[reihe[j]]*x[reihe[j+1]];
		j-=(D>2?D-1:0);

		if(j<1)return x[0]*x[0]*x[2];
		j-=1;

		if(j<1)return x[0]*x[2]*x[2];
		j-=1;

		if(j<1)return x[2]*x[2]*x[1];
		j-=1;

		if(j<1)return x[2]*x[1]*x[1];
		j-=1;

		if(j<1)
		{
			double product=1;
			for(int jj=0;jj<D;++jj)
				product*=x[reihe[jj]];
			return product;
		}
		j-=1;

		if (j < 1000) {
			double xx[D];
			for (int d = 0; d < D; ++d)
				xx[d] = (double) j / 1000. * 3. * x[d];
			return payoff(xx, zeit);
		}
		j -= 1000;

		if (j < 1000) {
			return max(y[reihe[0]]- (double)(j)/1000.*3,0);
		}
		j -= 1000;

		if (j < 1000) {
			double xx[D];
			for (int d = 0; d < D; ++d){
				xx[d] = x[d];
				if(j%D==d)
					xx[d] *=(double) j / 1000. * 6.;
			}
			return payoff(xx, zeit);
		}
		j -= 1000;

		if (j < 1000) {
			double xx[D];
			for (int d = 0; d < D; ++d){
				xx[d] = x[d];
				if(j%D!=d)
					xx[d] *=(double) j / 1000. * 6.;
			}
			return payoff(xx, zeit);
		}
		j -= 1000;

		if (j < 1000) {
			double diff = x[reihe[0]]-x[reihe[1]];
			return max(diff - (double) (j) / 1000. * 1. * (double) D * X0[0], 0);
		}
		j -= 1000;

		if (j < 1000) {
			double summe = 0;
			for (int d = 0; d < D; ++d)
				summe += x[d];
			return max(summe - (double) (j) / 1000. * 4. * (double) D * X0[0], 0);
		}
		j -= 1000;

		if (j < 1000) {
			double summe = 0;
			for (int d = 0; d < D; ++d)
				if(d!=j%3)
					summe += x[d];
			return max(summe - (double) (j) / 1000. * 3. * (double) D * X0[0], 0);
		}
		j -= 1000;
	}

	printf("Error598");
	return -1;
}
