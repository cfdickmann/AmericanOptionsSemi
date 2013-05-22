/*
 * Polynom.cpp
 *
 *  Created on: May 21, 2013
 *      Author: cfdickmann
 */

#include "Polynom.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

namespace std {

int binomial(int n, int k) {
	int num, den;
	if (n < k) {
		return (0);
	} else {
		den = 1;
		num = 1;
		for (int i = 1; i <= k; i = i + 1)
			den = den * i;
		for (int j = n - k + 1; j <= n; j = j + 1)
			num = num * j;
		return (num / den);
	}
}

Polynom potenzverschieben(int p, double v) {
	//(x+v)^k
	double erg[p + 1];
	for (int l = 0; l <= p; ++l)
		erg[l] = 0;

	for (int l = 0; l <= p; ++l) {
		erg[l] += (double) binomial(p, l) * pow(v, p - l);
//	printf("binomial %f, pow %f\n",(double)binomial(p, l),pow(v, p - l));
	}
//	return erg;

//	for (int l = 0; l <= p; ++l)
//		printf(",%f", erg[l]);
	Polynom e;
	e.set(erg, p + 1);
//	printf("hier:");
//	e.ausgeben();
//	delete[] erg;
	return e;
}

Polynom::Polynom() {
	length = 1;
	koeff = new double[1];
	koeff[0] = 0;
}

void Polynom::ausgeben() {
	printf("%.12lfx^%d\t", koeff[length - 1], length - 1);
	for (int k = length - 2; k >= 0; --k)
		printf(" + %.12lfx^%d\t", koeff[k], k);
	printf("\n");
}

void Polynom::malnehmen(double c) {
	for (int k = 0; k < length; ++k)
		koeff[k] *= c;
}

double Polynom::auswerten(double u) {
	double erg = 0;
	for (int k = 0; k < length; k++)
		erg += pow(u, k) * koeff[k];
	return erg;
}

void Polynom::multiply_with(Polynom* p) {
	double erg[length + p->length - 1];
	for (int k = 0; k < length + p->length - 1; ++k)
		erg[k] = 0;
	for (int kx = 0; kx < length; kx++)
		for (int ky = 0; ky < p->length; ++ky)
			erg[ky + kx] += koeff[kx] * p->koeff[ky];

	delete[] koeff;
	//koeff = erg;
	length += p->length - 1;
	koeff = new double[length];
	for (int ll = 0; ll < length; ++ll)
		koeff[ll] = erg[ll];
}

void Polynom::lerweitern(int l) {
	if (l <= length)
		return;

	double* kopie = new double[l];
	for (int ll = 0; ll < length; ++ll)
		kopie[ll] = koeff[ll];
	for (int ll = length; ll < l; ++ll)
		kopie[ll] = 0;

	delete[] koeff;
	koeff = new double[l];
	length = l;
	for (int ll = 0; ll < l; ++ll)
		koeff[ll] = kopie[ll];

	delete[] kopie;
}

void Polynom::add(Polynom* b) {
	if (b->length > length)
		lerweitern(b->length);
	for (int ll = 0; ll < length; ++ll)
		koeff[ll] += b->koeff[ll];
}

void Polynom::verschieben(double v) {
////		Polynom* qqq=potenzverschieben(1,2).ausgeben();
//	potenzverschieben(1, 2).ausgeben();
//	potenzverschieben(1, 2).ausgeben();
//	potenzverschieben(1, 2).ausgeben();
//	exit(0);

	Polynom e;
	double eins[1] = { koeff[0] };
	e.set(eins, 1);
//	e.ausgeben();
//e.lerweitern(10);
	for (int k = 1; k < length; ++k) {
//		printf("k=%d, v=%f\n",k,v);
		Polynom p = potenzverschieben(k, v);
//p.lerweitern(100);
		p.malnehmen(koeff[k]);
//		printf("Teilpolynom %d: ", k);
//		p.ausgeben();

//		p.ausgeben();
		e.add(&p);
//		delete p;
//		printf("Summe e : ");
//		e.ausgeben();
	}

	delete[] koeff;
	koeff = new double[e.length];
	length = e.length;
	for (int ll = 0; ll < e.length; ++ll)
		koeff[ll] = e.koeff[ll];

//	printf("\n");
//	Polynom p=potenzverschieben(,4);
////	p.lerweitern(7);
//	p.ausgeben();
//
//	Polynom q=potenzverschieben(5,1);
////	q.lerweitern(7);
//	q.ausgeben();
//
//	p.add(q);
//
//
//	printf("summe:\n");
//	p.ausgeben();
//
//	exit(0);
//printf("laenge: %d\n\n\n",p.length);

}

void Polynom::set(double *a, int l) {
	delete[] koeff;
	koeff = new double[l];
	length = l;
	for (int ll = 0; ll < l; ++ll)
		koeff[ll] = a[ll];
}

Polynom::~Polynom() {
	delete[] koeff;
}

} /* namespace std */
