/*
 *  Abeles.h
 *  Abeles
 *
 Copyright Andrew Nelson and ANSTO 2009
 *
 */
 //we can do the calculation multithreaded, it should be faster.
#ifdef MACIGOR
#include <pthread.h>
#endif
#ifdef WINIGOR
#include "pthread.h"
#include "sched.h"
#include "semaphore.h"
#endif

#include <complex>
#include "GaussWeights.h"

using namespace std;
// Prototypes
#ifdef MACIGOR
//using namespace std;
#include <sys/sysctl.h>
#endif	

extern GaussWeights *pinstance;
extern pthread_mutex_t changeWeightsMutex;
extern int NUM_CPUS;

HOST_IMPORT int XOPMain(IORecHandle ioRecHandle);
// Custom error codes
#define REQUIRES_IGOR_700 1 + FIRST_XOP_ERR
#define NON_EXISTENT_WAVE 2 + FIRST_XOP_ERR
#define REQUIRES_SP_OR_DP_WAVE 3 + FIRST_XOP_ERR
#define INCORRECT_INPUT 4 + FIRST_XOP_ERR
#define WAVES_NOT_SAME_LENGTH 5 + FIRST_XOP_ERR
#define REQUIRES_DP_WAVE 6 + FIRST_XOP_ERR

 int Abelescalc(const double*,double, double*);
 int Abelescalc_imag(const double*,double,double*);
 void matmul(complex<double>[2][2], complex<double>[2][2], complex<double>[2][2] );
 complex<double> fres(complex<double>, complex<double>, double);
