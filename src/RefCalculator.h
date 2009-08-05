/*
 *  Threadcalculator.h
 *  Abeles
 *
 *  Created by andrew on 24/02/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
	void *AbelesThreadWorker(void *arg);
	void *AbelesImagThreadWorker(void *arg);
//	using namespace std;



int AbelesCalcAll(double *coefP, double *yP, double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep);
int AbelesCalc_ImagAll(double *coefP, double *yP, double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep);

typedef struct{
	long npoints;
	int Vmullayers;
	int Vappendlayer;
	int Vmulrep;
	double *coefP;
	double *yP;
	double *xP;	
}  refCalcParm;

