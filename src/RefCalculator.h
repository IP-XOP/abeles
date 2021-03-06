/*
 *  Threadcalculator.h
 *  Abeles
 *
 Copyright Andrew Nelson and ANSTO 2009
 *
 */
	void *AbelesThreadWorker(void *arg);
	void *AbelesImagThreadWorker(void *arg);
//	using namespace std


int AbelesCalcAll(const double *coefP, double *yP, const double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep);
int AbelesCalc_ImagAll(const double *coefP, double *yP, const double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep);
int smearedAbelescalcAll(const double *coefP, double *yP, const double *xP, const double *dxP, long npoints);

//a structure that described the variables required to calculate the reflectivity
typedef struct{
	//number of Q points we have to calculate
	long npoints;
	//how many layers in the multilayer (optional, but if not specified should be 0)
	int Vmullayers;
	//where the multilayer is appended to the basic model.
	int Vappendlayer;
	//how many repeats in the multilayer
	int Vmulrep;
	//a double array containing the model coefficients (assumed to be correct)
	const double *coefP;
	//the Reflectivity values to return
	double *yP;
	//the Q values to do the calculation for.
	const double *xP;	
}  refCalcParm;

#define PI 3.14159265358979323846
