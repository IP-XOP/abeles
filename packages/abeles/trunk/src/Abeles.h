/*
 *  Abeles.h
 *  Abeles
 *
 *  Created by andrew on 24/02/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "MyComplex.h"

using namespace MyComplexNumber;
// Prototypes
#ifdef _MACINTOSH_
using namespace std;
#include <sys/sysctl.h>
HOST_IMPORT int main(IORecHandle ioRecHandle);
#endif	
#ifdef _WINDOWS
HOST_IMPORT void main(IORecHandle ioRecHandle);
#endif

// Custom error codes
#define REQUIRES_IGOR_400 1 + FIRST_XOP_ERR
#define NON_EXISTENT_WAVE 2 + FIRST_XOP_ERR
#define REQUIRES_SP_OR_DP_WAVE 3 + FIRST_XOP_ERR
#define INCORRECT_INPUT 4 + FIRST_XOP_ERR
#define WAVES_NOT_SAME_LENGTH 5 + FIRST_XOP_ERR
#define REQUIRES_DP_WAVE 6 + FIRST_XOP_ERR

#define PI 3.14159265358979323846

 int Abelescalc(double*,double, double*);
 int Abelescalc_imag(double*,double,double*);
 int smearedAbelescalcAll(double *coefP, double *yP, double *xP, double *dxP, long npoints);
 void matmul(MyComplex,MyComplex,MyComplex);
 MyComplex fres(MyComplex,MyComplex,double);