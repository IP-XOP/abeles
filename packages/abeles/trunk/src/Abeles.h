/*
 *  Abeles.h
 *  Abeles
 *
 Copyright Andrew Nelson and ANSTO 2009
 *
 */
#include "MyComplex.h"

using namespace MyComplexNumber;
// Prototypes
#ifdef MACIGOR
using namespace std;
#include <sys/sysctl.h>
#endif	

HOST_IMPORT int main(IORecHandle ioRecHandle);
// Custom error codes
#define REQUIRES_IGOR_600 1 + FIRST_XOP_ERR
#define NON_EXISTENT_WAVE 2 + FIRST_XOP_ERR
#define REQUIRES_SP_OR_DP_WAVE 3 + FIRST_XOP_ERR
#define INCORRECT_INPUT 4 + FIRST_XOP_ERR
#define WAVES_NOT_SAME_LENGTH 5 + FIRST_XOP_ERR
#define REQUIRES_DP_WAVE 6 + FIRST_XOP_ERR

 int Abelescalc(const double*,double, double*);
 int Abelescalc_imag(const double*,double,double*);
 void matmul(MyComplex,MyComplex,MyComplex);
 MyComplex fres(MyComplex,MyComplex,double);