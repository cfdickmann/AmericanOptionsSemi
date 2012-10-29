/*
 * NummerUndZeiger.h
 *
 *  Created on: Oct 29, 2012
 *      Author: cfdickmann
 */

#ifndef NUMMERUNDZEIGER_H_
#define NUMMERUNDZEIGER_H_

class ThreadDaten {
public:
	double* ergebnis;
	int nummer;
	void setErgebnis(double D);
	double getErgebnis();

	ThreadDaten();
	virtual ~ThreadDaten();
};

#endif /* NUMMERUNDZEIGER_H_ */
