/*	Abeles.cpp
calculates specular reflectivity as a function of Q momentum transfer.
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#include "Abeles.h"
#include <math.h>
#include <exception>
#include "RefCalculator.h"

//we can do the calculation multithreaded, it should be faster.
#ifdef MACIGOR
#include <pthread.h>
#endif
#ifdef WINIGOR
#include "pthread.h"
#include "sched.h"
#include "semaphore.h"
#endif

#pragma pack(2)	// All structures passed to Igor are two-byte aligned.

/*
this definition is for a fit function that takes a model wave and a Q value (passed as a variable) and returns
a single reflectivity value
*/
typedef struct FitParams {
	double x;				// Independent variable.
	waveHndl waveHandle;	// Coefficient wave.
	UserFunctionThreadInfoPtr tp;
	double result;
} FitParams, *FitParamsPtr;

/*
 this definition is for a fit function that takes a model wave and a Q wave (passed as a variable) and returns 
 reflectivity values, as a wave.
 */

typedef struct FitParamsAll {
	waveHndl XWaveHandle;	// X wave (input).
	waveHndl YWaveHandle;	// Y wave (output).
	waveHndl CoefHandle;	// Coefficient wave.
	UserFunctionThreadInfoPtr tp;
	double result;			// not actually used.
}FitParamsAll, *FitParamsAllPtr;

/*
 this definition is for a fit function that takes a model wave and a Q wave (passed as a variable) and 
 resolution in Q, dQ and returns 
 reflectivity values, as a wave.
*/

typedef struct SmearedParamsAll {
	waveHndl dXWaveHandle;	// dQ wave (input).	
	waveHndl XWaveHandle;	// Q wave (input).
	waveHndl YWaveHandle;	// R wave (output).
	waveHndl CoefHandle;	// Coefficient wave.
	UserFunctionThreadInfoPtr tp;
	double result;			// not actually used.
}SmearedParamsAll, *SmearedParamsAllPtr;

#pragma pack()	// All structures passed to Igor are two-byte aligned.


/*
 because we are going to do the calculation in a threaded fashion we need to know the number of CPU's.
 This is stored in a global so we don't need to find it out every time we call the function.
 */
int NUM_CPUS = 1;


extern "C" int
smearedAbelesAll(SmearedParamsAllPtr p){
	CountInt ncoefs, npoints;
	double realVal, imagVal;
	int nlayers, Vmullayers=-1, err=0;
	double *coefP = NULL;
	double *xP = NULL;
	double *yP = NULL;
	double *dxP = NULL;
	
	
	if (p->CoefHandle == NULL || p->YWaveHandle == NULL || p->XWaveHandle == NULL || p->dXWaveHandle == NULL) 
	{
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	if (WaveType(p->CoefHandle) != NT_FP64 ||
		  WaveType(p->YWaveHandle) != NT_FP64 ||
		   WaveType(p->XWaveHandle) != NT_FP64 ||
			WaveType(p->dXWaveHandle) != NT_FP64){
		SetNaN64(&p->result);
		err = REQUIRES_DP_WAVE;
		goto done;
	}
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1){
		Vmullayers=0;
	} else{
		Vmullayers=(long)realVal;
	}
	
	ncoefs= WavePoints(p->CoefHandle);
	npoints = WavePoints(p->YWaveHandle);
	if (npoints != WavePoints(p->XWaveHandle)){
		SetNaN64(&p->result);
		err = WAVES_NOT_SAME_LENGTH;
		goto done;
	}
	if (npoints != WavePoints(p->dXWaveHandle)){
		SetNaN64(&p->result);
		err = WAVES_NOT_SAME_LENGTH;
		goto done;
	}
	
	try{
		coefP =  new double[ncoefs];
		xP =  new double[npoints];
		yP = new double[npoints];
		dxP = new double[npoints];
	} catch (...){
		err = NOMEM;
		goto done;
	}
	
	if(err = MDGetDPDataFromNumericWave(p->CoefHandle, coefP))
		goto done;
	if(err = MDGetDPDataFromNumericWave(p->YWaveHandle, yP))
		goto done;
	if(err = MDGetDPDataFromNumericWave(p->XWaveHandle, xP))
		goto done;
	if(err = MDGetDPDataFromNumericWave(p->dXWaveHandle, dxP))
		goto done;
	
	nlayers = (long)coefP[0];
	if(ncoefs != (4*Vmullayers+(4*nlayers+6))){
		err = INCORRECT_INPUT;
		goto done;
	};
	
	if(err = smearedAbelescalcAll(coefP, yP, xP, dxP, (long) npoints))
		goto done;
	if(err = MDStoreDPDataInNumericWave(p->YWaveHandle, yP))
		goto done;
	
	WaveHandleModified(p->YWaveHandle);
	p->result = 0;		// not actually used by FuncFit
	
done:
	if(xP != NULL)
		delete [] xP;
	if(yP != NULL)
		delete [] yP;
	if(coefP != NULL)
		delete [] coefP;
	if(dxP != NULL)
		delete [] dxP;
	
	return err;	
}


void getMultiLayerParams(long *Vmullayers, long *Vappendlayer, long *Vmulrep){
	//temporary places for calculation
	double realVal, imagVal;
	
	//this may meanthere is a multilayer model specified
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1)
		*Vmullayers = 0;
	else
		*Vmullayers=(long)realVal;
	
	if(FetchNumVar("Vappendlayer", &realVal, &imagVal) == -1)
		*Vappendlayer = 0;
	else
		*Vappendlayer=(long)realVal;
	
	if(FetchNumVar("Vmulrep", &realVal, &imagVal) == -1)
		*Vmulrep=0;
	else
		*Vmulrep=(long)realVal;
	
}

extern "C"
int Abeles_bmagAll(FitParamsAllPtr p){
	int err = 0;
	
	//number of coefficients, number of Q points
	CountInt ncoefs,npoints;
	
	//how many layers you have
	long nlayers;
	
	//how many layers in multilayer (if you have one)
	long Vmullayers=-1;
	//where in the normal model the multilayer gets appended to
	long Vappendlayer=0;
	//how many times the multilayer repeats
	long Vmulrep=0;
	
	//a variable for iterating for loops
	long ii;
	long threadUsed;
	
	//we will extract values from the supplied waves and store them as double arrays
	//these pointers refer to the values extracted.
	double *coefP = NULL;
	double *xP = NULL;
	double *yP = NULL;

	double *tempYPplusplus = NULL;
	double *tempYPminusminus = NULL;
	double *coefPplusplus = NULL;
	double *coefPminusminus = NULL;
	
	
	//variables for the threadwise calculation of the reflectivity
	extern int NUM_CPUS;
	int threadsToCreate = 1;
	pthread_t *threads = NULL;
	refCalcParm *arg = NULL;
	long pointsEachThread = 0;
	long pointsRemaining = 0;
	long pointsConsumed = 0;
	
	//stuff starts here
	if (p->CoefHandle == NIL ||	p->YWaveHandle == NIL || p->XWaveHandle == NIL ){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	if (WaveType(p->CoefHandle) != NT_FP64 ||
		WaveType(p->YWaveHandle) != (NT_FP64 | NT_CMPLX) ||
		WaveType(p->XWaveHandle) != NT_FP64){
		SetNaN64(&p->result);
		err = REQUIRES_DP_WAVE;
		goto done;
	}
	
	ncoefs = WavePoints(p->CoefHandle);
	npoints = WavePoints(p->YWaveHandle);
	if (npoints != WavePoints(p->XWaveHandle)){
		SetNaN64(&p->result);
		err = WAVES_NOT_SAME_LENGTH;
		goto done;
	}
	
	coefP = (double*) WaveData(p->CoefHandle);
	yP = (double*) WaveData(p->YWaveHandle);
	xP = (double*) WaveData(p->XWaveHandle);
		
	//the number of layers should relate to the size of the coefficient wave (4*N+8 + 4*Vmullayers)
	//if not, the coefficient wave is wrong.
	nlayers = (long)coefP[0];
	
	if(ncoefs != 4 * nlayers + 8){
		//this may meanthere is a multilayer model specified
		getMultiLayerParams(&Vmullayers, &Vappendlayer, &Vmulrep);
				
		if(ncoefs != (4 * nlayers + 8) + 4 * Vmullayers){
			err = INCORRECT_INPUT;
			goto done;
		}				
	};
	
	try{
		coefPplusplus = new double[ncoefs - 2];
		coefPminusminus = new double[ncoefs -2];
		tempYPplusplus = new double[npoints];
		tempYPminusminus = new double[npoints];
	} catch (...){
		err = NOMEM;
		goto done;
	}
	
	coefPplusplus[0] = nlayers;
	coefPplusplus[1] = coefP[1];
	coefPplusplus[2] = coefP[2] + coefP[3];
	coefPplusplus[3] = coefP[4] + coefP[5];
	coefPplusplus[4] = coefP[6];
	coefPplusplus[5] = coefP[7];
	
	coefPminusminus[0] = nlayers;
	coefPminusminus[1] = coefP[1];
	coefPminusminus[2] = coefP[2] - coefP[3];
	coefPminusminus[3] = coefP[4] - coefP[5];
	coefPminusminus[4] = coefP[6];
	coefPminusminus[5] = coefP[7];
	
	for(ii = 0 ; ii < nlayers ; ii += 1){
		coefPplusplus[4 * ii + 6] = coefP[4 * ii + 8];
		coefPplusplus[4 * ii + 7] = coefP[4 * ii + 9] + coefP[4 * ii + 10];
		coefPplusplus[4 * ii + 8] = 0;
		coefPplusplus[4 * ii + 9] = coefP[4 * ii + 11];

		coefPminusminus[4 * ii + 6] = coefP[4 * ii + 8];
		coefPminusminus[4 * ii + 7] = coefP[4 * ii + 9] - coefP[4 * ii + 10];
		coefPminusminus[4 * ii + 8] = 0;
		coefPminusminus[4 * ii + 9] = coefP[4 * ii + 11];
	}
	
	threadsToCreate = NUM_CPUS;
	
	//create threads for the calculation
	threads = (pthread_t *) malloc((threadsToCreate) * sizeof(pthread_t));
	if(!threads){
		err = NOMEM;
		goto done;
	}
	
	//create arguments to be supplied to each of the threads
	arg = (refCalcParm *) malloc (sizeof(refCalcParm) * (threadsToCreate));
	if(!arg){
		err = NOMEM;
		goto done;
	}
	
	//need to calculated how many points are given to each thread.
	pointsEachThread = floorl(npoints / (threadsToCreate / 2L));
	pointsRemaining = npoints;
	pointsConsumed = 0;
	threadUsed = 0;
	
	for (ii = 0; ii < (long)(threadsToCreate / 2L); ii++){
		if(ii == (long)(threadsToCreate / 2L) - 1)
		   pointsEachThread = pointsRemaining;
		   
		arg[threadUsed].coefP = coefPplusplus;
		arg[threadUsed].npoints = (long) pointsEachThread;
		
		arg[threadUsed].Vmullayers = Vmullayers;
		arg[threadUsed].Vappendlayer = Vappendlayer;
		arg[threadUsed].Vmulrep = Vmulrep;
		
		//the following two lines specify where the Q values and R values will be sourced/written.
		//i.e. an offset of the original array.
		arg[threadUsed].xP = xP + pointsConsumed;
		arg[threadUsed].yP = tempYPplusplus + pointsConsumed;
		
		pthread_create(&threads[threadUsed], NULL, AbelesThreadWorker, (void *)(arg + threadUsed));
		pointsRemaining -= pointsEachThread;
		pointsConsumed += pointsEachThread;
		threadUsed += 1;
	}

	//need to calculated how many points are given to each thread.
   pointsEachThread = floorl(npoints / (threadsToCreate / 2L));
   pointsRemaining = npoints;
   pointsConsumed = 0;
   
   for (ii = 0; ii < (long)(threadsToCreate / 2L); ii++){
	   if(ii == (long)(threadsToCreate / 2L) - 1)
		  pointsEachThread = pointsRemaining;
		  
		  arg[threadUsed].coefP = coefPminusminus;
		  arg[threadUsed].npoints = (long) pointsEachThread;
		  
		  arg[threadUsed].Vmullayers = Vmullayers;
		  arg[threadUsed].Vappendlayer = Vappendlayer;
		  arg[threadUsed].Vmulrep = Vmulrep;
		  
		  //the following two lines specify where the Q values and R values will be sourced/written.
		  //i.e. an offset of the original array.
		  arg[threadUsed].xP = xP + pointsConsumed;
		  arg[threadUsed].yP = tempYPminusminus + pointsConsumed;
		  
		  pthread_create(&threads[threadUsed], NULL, AbelesThreadWorker, (void *)(arg + threadUsed));
		  pointsRemaining -= pointsEachThread;
		  pointsConsumed += pointsEachThread;
		  threadUsed += 1;
	}
				  
	for (ii = 0; ii < threadsToCreate; ii++)
		pthread_join(threads[ii], NULL);

	for(ii = 0 ; ii < npoints ; ii += 1){
		yP[2 * ii] = tempYPplusplus[ii]; 
		yP[2 * ii + 1] = tempYPminusminus[ii]; 
	}

	WaveHandleModified(p->YWaveHandle);
	p->result = 0;		// not actually used by FuncFit
	
done:
	if(tempYPplusplus)
		delete [] tempYPplusplus;
	if(tempYPminusminus)
		delete [] tempYPminusminus;
	if(coefPminusminus)
		delete [] coefPminusminus;
	if(coefPplusplus)
		delete [] coefPplusplus;
	
	if(threads)
		free(threads);
	if(arg)
		free(arg);
	
	return err;			
}

/*mode:
0 = AbelesAll
1 = Abeles_imagAll
2 = Abeles_BmagAll
 */
 
int AbelesAllWrapper(FitParamsAllPtr p, int mode){
	int err = 0;

	//number of coefficients, number of Q points
	CountInt ncoefs,npoints;
	
	//how many layers you have
	long nlayers;
	
	//how many layers in multilayer (if you have one)
	long Vmullayers=-1;
	//where in the normal model the multilayer gets appended to
	long Vappendlayer=0;
	//how many times the multilayer repeats
	long Vmulrep=0;
	
	//a variable for iterating for loops
	long ii;
	
	//we will extract values from the supplied waves and store them as double arrays
	//these pointers refer to the values extracted.
	double *coefP = NULL;
	double *xP = NULL;
	double *yP = NULL;
	
	//variables for the threadwise calculation of the reflectivity
	extern int NUM_CPUS;
	int threadsToCreate = 1;
	double isItWorthThreading = 0;
	pthread_t *threads = NULL;
	refCalcParm *arg = NULL;
	long pointsEachThread = 0;
	long pointsRemaining = 0;
	long pointsConsumed = 0;
	
	//the functions that will do the calculation
	void* (*threadWorkerFunc)(void*);
	int (*calcAllFunc)(const double *, double *, const double *,long, int, int, int);
	
	//this is now threadsafe, check IGOR is >6.20
	if (igorVersion < 620)
		if (!RunningInMainThread())
			return NOT_IN_THREADSAFE;
		
	//stuff starts here
	if (p->CoefHandle == NIL ||	p->YWaveHandle == NIL || p->XWaveHandle == NIL ){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	if (WaveType(p->CoefHandle) != NT_FP64 ||
		WaveType(p->YWaveHandle) != NT_FP64 ||
		WaveType(p->XWaveHandle) != NT_FP64){
		SetNaN64(&p->result);
		err = REQUIRES_DP_WAVE;
		goto done;
	}
	
	ncoefs = WavePoints(p->CoefHandle);
	npoints = WavePoints(p->YWaveHandle);
	if (npoints != WavePoints(p->XWaveHandle)){
		SetNaN64(&p->result);
		err = WAVES_NOT_SAME_LENGTH;
		goto done;
	}
	
	coefP = (double*) WaveData(p->CoefHandle);
	yP = (double*) WaveData(p->YWaveHandle);
	xP = (double*) WaveData(p->XWaveHandle);
	
	//the number of layers should relate to the size of the coefficient wave (4*N+6 + 4*Vmullayers)
	//if not, the coefficient wave is wrong.
	nlayers = (long)coefP[0];
	
	switch (mode) {
		case 0:
			threadWorkerFunc = &AbelesThreadWorker;
			calcAllFunc = &AbelesCalcAll;

			if(ncoefs != 4 * nlayers + 6){				
				//this may meanthere is a multilayer model specified
				getMultiLayerParams(&Vmullayers, &Vappendlayer, &Vmulrep);
				
				if(ncoefs != (4 * nlayers + 6) + 4 * Vmullayers){
					err = INCORRECT_INPUT;
					goto done;
				}
				
			};
			
			break;
		case 1:
			threadWorkerFunc = &AbelesImagThreadWorker;
			calcAllFunc = &AbelesCalc_ImagAll;
			
			if(ncoefs != 4 * nlayers + 8){
				//this may meanthere is a multilayer model specified
				getMultiLayerParams(&Vmullayers, &Vappendlayer, &Vmulrep);
				
				if(ncoefs != (4 * nlayers + 8) + 4 * Vmullayers){
					err = INCORRECT_INPUT;
					goto done;
				}				
			};
			
			break;
		default:
			break;
	}
		
	//this relationship was worked out for a dualcore machine.
	//I worked out how long it took for a certain number of points and a certain number of layers
	//then I calculated the cross over for when it was worth threading, rather than not.
	//i.e. for a given number of layers, how many points were required for multithreading to be worthwhile.
	//I plotted the number of points (y) vs the number of layers (x), giving the following relationship.
	isItWorthThreading = 3.382 + 641. *pow(coefP[0], -0.73547);
	if((double) npoints < isItWorthThreading)
		threadsToCreate = 1;
	else
		threadsToCreate = NUM_CPUS;
	
	//create threads for the calculation
	threads = (pthread_t *) malloc((threadsToCreate - 1) * sizeof(pthread_t));
	if(!threads && NUM_CPUS > 1){
		err = NOMEM;
		goto done;
	}
	
	//create arguments to be supplied to each of the threads
	arg = (refCalcParm *) malloc (sizeof(refCalcParm)*(threadsToCreate - 1));
	if(!arg && NUM_CPUS > 1){
		err = NOMEM;
		goto done;
	}
	
	//need to calculated how many points are given to each thread.
	pointsEachThread = floorl(npoints / threadsToCreate);
	pointsRemaining = npoints;
	
	for (ii = 0; ii < threadsToCreate - 1; ii++){
		arg[ii].coefP = coefP;
		arg[ii].npoints = pointsEachThread;
		
		arg[ii].Vmullayers = Vmullayers;
		arg[ii].Vappendlayer = Vappendlayer;
		arg[ii].Vmulrep = Vmulrep;
		
		//the following two lines specify where the Q values and R values will be sourced/written.
		//i.e. an offset of the original array.
		arg[ii].xP = xP + pointsConsumed;
		arg[ii].yP = yP + pointsConsumed;
		
		pthread_create(&threads[ii], NULL, threadWorkerFunc, (void *)(arg+ii));
		pointsRemaining -= pointsEachThread;
		pointsConsumed += pointsEachThread;
	}
	
	//do the remaining points in the main thread.
	if(err = calcAllFunc(coefP, yP+pointsConsumed, xP+pointsConsumed, pointsRemaining, Vmullayers, Vappendlayer, Vmulrep))
		goto done;
	
	for (ii = 0; ii < threadsToCreate - 1 ; ii++)
		pthread_join(threads[ii], NULL);
	
	WaveHandleModified(p->YWaveHandle);
	p->result = 0;		// not actually used by FuncFit
	
done:
	if(threads)
		free(threads);
	if(arg)
		free(arg);
	
	return err;		
		
}

/*
 An all at once fit function.  You send all the Qvalues and a model wave, and get all the reflectivity values back.
 */
extern "C" int
AbelesAll(FitParamsAllPtr p){
	return AbelesAllWrapper(p, 0);
}

extern "C" int
Abeles_imagAll(FitParamsAllPtr p){
	return AbelesAllWrapper(p, 1);
}

extern "C" int
Abeles(FitParamsPtr p){
	int np, err = 0;
	double *Abelesparams = NULL;
	double x;
	int Vmullayers = 0;

	//this is now threadsafe, check IGOR is >6.20
	if (igorVersion < 620)
		if (!RunningInMainThread())
			return NOT_IN_THREADSAFE;
	
	if (p->waveHandle == NULL){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	
	if(WaveType(p->waveHandle) != NT_FP64){
		err = REQUIRES_DP_WAVE;
		SetNaN64(&p->result);
		goto done;
	}
			
	np= WavePoints(p->waveHandle);
	
	x = p->x;
	Abelesparams = (double*) WaveData(p->waveHandle);

	Vmullayers = ((np - 6) / 4) - (int) Abelesparams[0];
	
	if((int)np!=(int)(4 * Vmullayers + 4 * Abelesparams[0] + 6)){
		err = INCORRECT_INPUT;
		goto done;
	};
	
	if(err = Abelescalc(Abelesparams, x, &p->result))
		goto done;
done:
	
	return err;
}

extern "C" int
Abeles_imag(FitParamsPtr p){
	int np, err = 0;
	double *Abelesparams = NULL;
	double x,result;
	int Vmullayers = 0;
	
	//this is now threadsafe, check IGOR is >6.20
	if (igorVersion < 620)
		if (!RunningInMainThread())
			return NOT_IN_THREADSAFE;
	
	if (p->waveHandle == NULL){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	
	if(WaveType(p->waveHandle) != NT_FP64){
		err = REQUIRES_DP_WAVE;
		SetNaN64(&p->result);
		goto done;
	}
	
	np= WavePoints(p->waveHandle);
	Abelesparams= (double*) WaveData(p->waveHandle);
	x= p->x;
	
	Vmullayers = ((np - 8) / 4) - (int) Abelesparams[0];

	if((int)np != (int)(4 * Vmullayers + 4 * Abelesparams[0] + 8)){
		err = INCORRECT_INPUT;
		goto done;
	};
	
	if(err = Abelescalc_imag(Abelesparams, x, &result))
		goto done;
	p->result = result;
	
done:
	
	return err;
}

extern "C" int
parrattReflectance(FitParamsAllPtr p){
	//number of coefficients, number of Q points
	long ncoefs,npoints;
	
	//temporary places for calculation
	double realVal,imagVal;
	
	//how many layers you have
	long nlayers;
	
	//how many layers in multilayer (if you have one)
	long Vmullayers=-1;
	//where in the normal model the multilayer gets appended to
	long Vappendlayer=0;
	//how many times the multilayer repeats
	long Vmulrep=0;
	
	//the return code of the function
	long err=0;
	//a variable for iterating for loops
	long ii;
	
	//we will extract values from the supplied waves and store them as double arrays
	//these pointers refer to the values extracted.
	double *coefP = NULL;
	double *xP = NULL;
	double *yP = NULL;
	
	//variables for the threadwise calculation of the reflectivity
	extern int NUM_CPUS;
	int threadsToCreate = 1;
	double isItWorthThreading = 0;
	pthread_t *threads = NULL;
	refCalcParm *arg = NULL;
	long pointsEachThread = 0;
	long pointsRemaining = 0;
	long pointsConsumed = 0;
	
	//have to check that all the supplied wave references exist
	if (p->CoefHandle == NULL ||
		p->YWaveHandle == NULL ||
		p->XWaveHandle == NULL ){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	//check that all the supplied waves have the correct numerical type.
	if (!(WaveType(p->CoefHandle) != NT_FP64 ||
		   WaveType(p->YWaveHandle) != NT_FP64 ||
		    WaveType(p->XWaveHandle) != NT_FP64)){
		SetNaN64(&p->result);
		err = REQUIRES_DP_WAVE;
		goto done;
	}
	
	ncoefs= WavePoints(p->CoefHandle);
	npoints = WavePoints(p->YWaveHandle);
	if (npoints != WavePoints(p->XWaveHandle)){
		SetNaN64(&p->result);
		err = WAVES_NOT_SAME_LENGTH;
		goto done;
	}
	
	try{
		coefP =  new double[ncoefs];
		xP =  new double[npoints];
		yP = new double[npoints];
	} catch (...){
		err = NOMEM;
		goto done;
	}
	
	if(err = MDGetDPDataFromNumericWave(p->CoefHandle, coefP))
		goto done;
	if(err = MDGetDPDataFromNumericWave(p->YWaveHandle, yP))
		goto done;
	if(err = MDGetDPDataFromNumericWave(p->XWaveHandle, xP))
		goto done;
	
	//the number of layers should relate to the size of the coefficient wave (4*N+6 + 4*Vmullayers)
	//if not, the coefficient wave is wrong.
	nlayers = (long)coefP[0];
	if(ncoefs != 4 * nlayers + 6){
		//this may meanthere is a multilayer model specified
		if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1)
			Vmullayers=0;
		else
			Vmullayers=(long)realVal;
		
		if(ncoefs != (4 * nlayers + 6) + 4 * Vmullayers){
			err = INCORRECT_INPUT;
			goto done;
		}
		
		if(FetchNumVar("Vappendlayer", &realVal, &imagVal) == -1)
			Vappendlayer = 0;
		else
			Vappendlayer=(long)realVal;
		
		if(FetchNumVar("Vmulrep", &realVal, &imagVal) == -1)
			Vmulrep=0;
		else
			Vmulrep=(long)realVal;
	};
	
	//this relationship was worked out for a dualcore machine.
	//I worked out how long it took for a certain number of points and a certain number of layers
	//then I calculated the cross over for when it was worth threading, rather than not.
	//i.e. for a given number of layers, how many points were required for multithreading to be worthwhile.
	//I plotted the number of points (y) vs the number of layers (x), giving the following relationship.
	isItWorthThreading = 3.382 + 641. *pow(coefP[0], -0.73547);
	if((double) npoints < isItWorthThreading)
		threadsToCreate = 1;
	else
		threadsToCreate = NUM_CPUS;
	
	//create threads for the calculation
	threads = (pthread_t *) malloc((threadsToCreate - 1) * sizeof(pthread_t));
	if(!threads && NUM_CPUS > 1){
		err = NOMEM;
		goto done;
	}
	//create arguments to be supplied to each of the threads
	arg = (refCalcParm *) malloc (sizeof(refCalcParm)*(threadsToCreate - 1));
	if(!arg && NUM_CPUS > 1){
		err = NOMEM;
		goto done;
	}
	//need to calculated how many points are given to each thread.
	pointsEachThread = floorl(npoints / threadsToCreate);
	pointsRemaining = npoints;
	
	//if you have two CPU's, only create one extra thread because the main thread does half the work
	for (ii = 0; ii < threadsToCreate - 1; ii++){
		arg[ii].coefP = coefP;
		arg[ii].npoints = pointsEachThread;
		
		arg[ii].Vmullayers = Vmullayers;
		arg[ii].Vappendlayer = Vappendlayer;
		arg[ii].Vmulrep = Vmulrep;
		
		//the following two lines specify where the Q values and R values will be sourced/written.
		//i.e. an offset of the original array.
		arg[ii].xP = xP+pointsConsumed;
		arg[ii].yP = yP+pointsConsumed;
		
		pthread_create(&threads[ii], NULL, realReflectanceThreadWorker, (void *)(arg+ii));
		pointsRemaining -= pointsEachThread;
		pointsConsumed += pointsEachThread;
	}
	
	//do the remaining points in the main thread.
	if(err = realReflectance(coefP, yP+pointsConsumed, xP+pointsConsumed, pointsRemaining))
		goto done;
	
	
	for (ii = 0; ii < threadsToCreate - 1 ; ii++)
		pthread_join(threads[ii], NULL);
	
	//put the reflectivity data back into the y wave supplied.
	if(err = MDStoreDPDataInNumericWave(p->YWaveHandle, yP))
		goto done;
	
	WaveHandleModified(p->YWaveHandle);
	p->result = 0;		// not actually used by FuncFit
	
done:
	//now have to free the thread memory and argument memory
	if(threads)
		free(threads);
	if(arg)
		free(arg);
	
	if(xP != NULL)
		delete [] xP;
	if(yP != NULL)
		delete [] yP;
	if(coefP != NULL)
		delete [] coefP;
	
	return err;	
}


static XOPIORecResult
RegisterFunction()
{
	int funcIndex;
	
	funcIndex = GetXOPItem(0);			// Which function invoked ?
	switch (funcIndex) {
		case 0:							// y = Abeles(w,x) (curve fitting function).
			return((XOPIORecResult)Abeles);	// This function is called using the direct method.
			break;
		case 1:
			return((XOPIORecResult)Abeles_imag);
			break;
		case 2:
			return((XOPIORecResult)AbelesAll);
			break;
		case 3:
			return((XOPIORecResult)Abeles_imagAll);
			break;
		case 4:
			return((XOPIORecResult)smearedAbelesAll);
			break;
		case 5:
			return((XOPIORecResult)parrattReflectance);
			break;
		case 6:
			return((XOPIORecResult)Abeles_bmagAll);
	}
	return 0;
}

/*	XOPEntry()
 
 This is the entry point from the host application to the XOP for all
 messages after the INIT message.
 */
extern "C" void
XOPEntry(void)
{	
	XOPIORecResult result = 0;
	
	switch (GetXOPMessage()) {
		case FUNCADDRS:
			result = RegisterFunction();	// This tells Igor the address of our function.
			break;
	}
	SetXOPResult(result);
}

/*	main(ioRecHandle)
 
 This is the initial entry point at which the host application calls XOP.
 The message sent by the host must be INIT.
 main() does any necessary initialization and then sets the XOPEntry field of the
 ioRecHandle to the address to be called for future messages.
 */
HOST_IMPORT int main(IORecHandle ioRecHandle){	
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	extern int NUM_CPUS;

	//find out the number of CPU's.
#ifdef WINIGOR
    SYSTEM_INFO sysInfo;  
     GetSystemInfo(&sysInfo);  
    NUM_CPUS = sysInfo.dwNumberOfProcessors;  
#endif
	
#ifdef MACIGOR
	int mib[2];
	size_t len;
	
	mib[0] = CTL_HW;
	mib[1] = HW_NCPU;
	len = sizeof(NUM_CPUS);
	sysctl(mib, 2, &NUM_CPUS, &len, NULL, 0);
#endif
	
	if (igorVersion < 600){
		SetXOPResult(IGOR_OBSOLETE);
		return EXIT_FAILURE;
	}else{
		SetXOPResult(0L);
		return EXIT_SUCCESS;
	}
}


MyComplex fres(MyComplex a, MyComplex b, double rough){
	return (compexp(-2 * rough * rough * a * b)) * (a - b)/(a + b);
}

int
Abelescalc(const double *coefP, double x, double *result){
	int err = 0;
	
	int Vmulrep = 0,Vmulappend = 0, Vmullayers=0;
	double realVal, imagVal;
	register int ii=0, jj=0, kk=0;
	
	double scale,bkg,subrough;
	double num=0,den=0, answer=0,qq;
	double anum,anum2;
	MyComplex temp,SLD,beta,rj;
	double numtemp=0;
	int offset=0;
	MyComplex  MRtotal[2][2];
	MyComplex subtotal[2][2];
	MyComplex MI[2][2];
	MyComplex temp2[2][2];
	MyComplex qq2;
	MyComplex oneC = MyComplex(1,0);
	MyComplex *pj_mul = NULL;
	MyComplex *pj = NULL;
	double *SLDmatrix = NULL;
	double *SLDmatrixREP = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		pj = new MyComplex [nlayers+2];
		SLDmatrix = new double [nlayers+2];
	} catch(...){
		err = NOMEM;
		goto done;
	}
	
	memset(pj, 0, sizeof(pj));
	memset(SLDmatrix, 0, sizeof(SLDmatrix));
	
	scale = coefP[1];
	bkg = fabs(coefP[4]);
	subrough = coefP[5];
	
	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 6;
	
	//fillout all the SLD's for all the layers
	for(ii=1; ii<nlayers+1;ii+=1){
		numtemp = 1e-6 * ((100 - coefP[4*ii+4])/100) * coefP[4*ii+3]+ (coefP[4*ii+4]*coefP[3]*1e-6)/100;		//sld of the layer
		*(SLDmatrix+ii) = 4 * PI * (numtemp  - (coefP[2]*1e-6));
	}
	*(SLDmatrix) = 0;
	*(SLDmatrix+nlayers+1) = 4* PI * ((coefP[3]*1e-6) - (coefP[2]*1e-6));
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) != -1){ // Fetch value
		Vmullayers = (int)realVal;
		if(FetchNumVar("Vappendlayer", &realVal, &imagVal) != -1)
			Vmulappend = (int)realVal;
		if(FetchNumVar("Vmulrep", &realVal, &imagVal) != -1) // Fetch value
			Vmulrep = (int)realVal;
		
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
			//set up an array for wavevectors
			try{
				SLDmatrixREP = new double [Vmullayers];
				pj_mul = new MyComplex [Vmullayers];
			} catch(...){
				err = NOMEM;
				goto done;
			}
			memset(pj_mul, 0, sizeof(pj_mul));
			for(ii = 0 ; ii < Vmullayers ; ii += 1){
				numtemp = (coefP[3]*1e-6*coefP[(4*ii)+offset+2]/100) +(1e-6 * ((100 - coefP[(4*ii)+offset+2])/100) * coefP[(4*ii)+offset+1]);		//sld of the layer
				*(SLDmatrixREP + ii) = 4 * PI * (numtemp  - (coefP[2] * 1e-6));
			}
		}
	}
	
	//intialise the matrices
	memset(MRtotal, 0, sizeof(MRtotal));
	MRtotal[0][0] = oneC ; MRtotal[1][1] = oneC;
	
	qq = x*x/4;
	qq2=MyComplex(qq, 0);
	
	for(ii=0; ii<nlayers+2 ; ii++){			//work out the wavevector in each of the layers
		pj[ii] = (*(SLDmatrix+ii)>qq) ? compsqrt(qq2-MyComplex(*(SLDmatrix+ii),0)): MyComplex(sqrt(qq-*(SLDmatrix+ii)),0);
		//pj[ii] = (*(SLDmatrix+ii)>qq) ? oneC.compsqrt(qq2-MyComplex(*(SLDmatrix+ii),0)): MyComplex(sqrt(qq-*(SLDmatrix+ii)),0);
	}
	
	//workout the wavevector in the toplayer of the multilayer, if it exists.
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
		memset(subtotal,0,sizeof(subtotal));
		subtotal[0][0]=MyComplex(1,0);subtotal[1][1]=MyComplex(1,0);
		pj_mul[0] = (*(SLDmatrixREP)>qq) ? compsqrt(qq2-MyComplex(*SLDmatrixREP,0)): MyComplex(sqrt(qq-*SLDmatrixREP),0);
		//	pj_mul[0] = (*(SLDmatrixREP)>qq) ? oneC.compsqrt(qq2-MyComplex(*SLDmatrixREP,0)): MyComplex(sqrt(qq-*SLDmatrixREP),0);
	}
	
	//now calculate reflectivities
	for(ii = 0 ; ii < nlayers+1 ; ii++){
		//work out the fresnel coefficients
		//this looks more complicated than it really is.
		//the reason it looks so convoluted is because if there is no complex part of the wavevector,
		//then it is faster to do the calc with real arithmetic then put it into a complex number.
		if(Vmullayers>0 && ii==Vmulappend && Vmulrep>0 ){
			rj=fres(pj[ii],pj_mul[0],coefP[offset+3]);
		} else {
			if((pj[ii]).im == 0 && (pj[ii+1]).im==0){
				anum = (pj[ii]).re;
				anum2 = (pj[ii+1]).re;
				rj = (ii==nlayers) ?
				MyComplex(((anum-anum2)/(anum+anum2))*exp(anum*anum2*-2*subrough*subrough),0)
				:
				MyComplex(((anum-anum2)/(anum+anum2))*exp(anum*anum2*-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5]),0);
			} else {
				rj = (ii == nlayers) ?
				((pj[ii]-pj[ii+1])/(pj[ii]+pj[ii+1]))*compexp(pj[ii]*pj[ii+1]*MyComplex(-2*subrough*subrough,0))
				:
				((pj[ii]-pj[ii+1])/(pj[ii]+pj[ii+1]))*compexp(pj[ii]*pj[ii+1]*MyComplex(-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5],0));	
			};
		}
		
		//work out the beta for the (non-multi)layer
		beta = (ii==0)? oneC : compexp(pj[ii] * MyComplex(0,fabs(coefP[4*ii+2])));
		
		//this is the characteristic matrix of a layer
		MI[0][0]=beta;
		MI[0][1]=rj*beta;
		MI[1][1]=oneC/beta;
		MI[1][0]=rj*MI[1][1];
		
		temp2[0][0] = MRtotal[0][0];
		temp2[0][1] = MRtotal[0][1];
		temp2[1][0] = MRtotal[1][0];
		temp2[1][1] = MRtotal[1][1];
		//multiply MR,MI to get the updated total matrix.			
		matmul(temp2,MI,MRtotal);
		
		if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0){
			//workout the wavevectors in each of the layers
			for(jj=1 ; jj < Vmullayers; jj++){
				pj_mul[jj] = (*(SLDmatrixREP+jj)>qq) ? compsqrt(qq2-MyComplex(*(SLDmatrixREP+jj),0)): MyComplex(sqrt(qq-*(SLDmatrixREP+jj)),0);
			}
			
			//work out the fresnel coefficients
			for(jj = 0 ; jj < Vmullayers; jj++){
				
				rj = (jj == Vmullayers-1) ?
				//if you're in the last layer then the roughness is the roughness of the top
				((pj_mul[jj]-pj_mul[0])/(pj_mul[jj]+pj_mul[0]))*compexp((pj_mul[jj]*pj_mul[0])*MyComplex(-2*coefP[offset+3]*coefP[offset+3],0))
				:
				//otherwise it's the roughness of the layer below
				((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*MyComplex(-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3],0));
				
				//Beta's
				beta = compexp(MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]);
				
				MI[0][0]=beta;
				MI[0][1]=rj*beta;
				MI[1][1]=oneC/beta;
				MI[1][0]=rj*MI[1][1];
				
				temp2[0][0] = subtotal[0][0];
				temp2[0][1] = subtotal[0][1];
				temp2[1][0] = subtotal[1][0];
				temp2[1][1] = subtotal[1][1];
				
				matmul(temp2,MI,subtotal);
			};
			
			for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
				if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
					for(jj=0;jj<Vmullayers;jj++){
						beta = compexp((MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]));
						
						if(jj==Vmullayers-1){
							if(Vmulappend==nlayers){
								rj = ((pj_mul[Vmullayers-1]-pj[nlayers+1])/(pj_mul[Vmullayers-1]+pj[nlayers+1]))*compexp((pj_mul[Vmullayers-1]*pj[nlayers+1])*MyComplex(-2*subrough*subrough,0));
							} else {
								rj = ((pj_mul[Vmullayers-1]-pj[Vmulappend+1])/(pj_mul[Vmullayers-1]+pj[Vmulappend+1]))*compexp((pj_mul[Vmullayers-1]*pj[Vmulappend+1])*MyComplex(-2*coefP[4*(Vmulappend+1)+5]*coefP[4*(Vmulappend+1)+5],0));
							};
						} else {
							rj = ((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*MyComplex(-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3],0));
						}
						
						MI[0][0]=beta;
						MI[0][1]=(rj*beta);
						MI[1][1]=oneC/beta;
						MI[1][0]=(rj*MI[1][1]);
						
						temp2[0][0] = MRtotal[0][0];
						temp2[0][1] = MRtotal[0][1];
						temp2[1][0] = MRtotal[1][0];
						temp2[1][1] = MRtotal[1][1];
						
						matmul(temp2,MI,MRtotal);
					}
				} else {
					temp2[0][0] = MRtotal[0][0];
					temp2[0][1] = MRtotal[0][1];
					temp2[1][0] = MRtotal[1][0];
					temp2[1][1] = MRtotal[1][1];
					
					matmul(temp2,subtotal,MRtotal);
				};
			};
		};
		
	}
	
	den=compnorm(MRtotal[0][0]);
	num=compnorm(MRtotal[1][0]);
	answer=(num)/(den);
	answer=(answer*scale)+bkg;
	
	*result = answer;
	
done:
	if(pj != NULL)
		delete [] pj;
	if(pj_mul !=NULL)
		delete[] pj_mul;
	if(SLDmatrix != NULL)
		delete[] SLDmatrix;
	if(SLDmatrixREP != NULL)
		delete[] SLDmatrixREP;
	
	return err;
}


int
Abelescalc_imag(const double *coefP, double x, double *result){
	int err = 0;
	
	int Vmulrep=0,Vmulappend=0,Vmullayers=0;
	double realVal,imagVal;
	int ii=0,jj=0,kk=0;
	
	double scale,bkg,subrough;
	double num=0,den=0, answer=0;
	
	MyComplex super;
	MyComplex sub;
	MyComplex temp,SLD,beta,rj,arg;
	MyComplex oneC = MyComplex(1,0);
	int offset=0;
	MyComplex MRtotal[2][2];
	MyComplex subtotal[2][2];
	MyComplex MI[2][2];
	MyComplex temp2[2][2];
	MyComplex qq2;
	MyComplex *pj_mul = NULL;
	MyComplex *pj = NULL;
	MyComplex *SLDmatrix = NULL;
	MyComplex *SLDmatrixREP = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		pj = new MyComplex[nlayers+2];
		SLDmatrix = new MyComplex[nlayers+2];
	} catch(...){
		err = NOMEM;
		goto done;
	}
	
	memset(pj, 0, sizeof(pj));
	memset(SLDmatrix, 0, sizeof(SLDmatrix));
	
	scale = coefP[1];
	bkg = coefP[6];
	subrough = coefP[7];
	sub= MyComplex(coefP[4]*1e-6,coefP[5]);
	super = MyComplex(coefP[2]*1e-6,coefP[3]);
	
	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 8;
	
	//fillout all the SLD's for all the layers
	for(ii=1; ii<nlayers+1;ii+=1){
		*(SLDmatrix+ii) = MyComplex(4 * PI , 0)*(MyComplex(coefP[4*ii+5]*1e-6,coefP[4*ii+6])-super);
	}
	*(SLDmatrix) = MyComplex(0,0);
	*(SLDmatrix+nlayers+1) = MyComplex(4 * PI, 0) * (sub - super);
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal)!=-1){ // Fetch value
		Vmullayers=(int)realVal;
		if(FetchNumVar("Vappendlayer", &realVal, &imagVal)!=-1) // Fetch value
			Vmulappend=(int)realVal;
		if(FetchNumVar("Vmulrep", &realVal, &imagVal) !=-1) // Fetch value
			Vmulrep=(int)realVal;
		
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
			//set up an array for wavevectors
			try{
				SLDmatrixREP = new MyComplex[Vmullayers];
				pj_mul = new MyComplex[Vmullayers];
			} catch(...){
				err = NOMEM;
				goto done;
			}
			memset(pj_mul, 0, sizeof(pj_mul));
			memset(SLDmatrixREP,0,sizeof(SLDmatrixREP));
			for(ii=0; ii<Vmullayers;ii+=1){
				*(SLDmatrixREP+ii) = MyComplex(4 * PI, 0)*(MyComplex(coefP[(4*ii)+offset+1]*1e-6,coefP[(4*ii)+offset+2])  - super);
			}
		}
	}
	
	//intialise the matrices
	memset(MRtotal,0,sizeof(MRtotal));
	MRtotal[0][0]=oneC;MRtotal[1][1]=oneC;
	
	qq2=MyComplex(x*x/4,0);
	
	for(ii=0; ii<nlayers+2 ; ii++){			//work out the wavevector in each of the layers
		pj[ii] = compsqrt(qq2-*(SLDmatrix+ii));
	}
	
	//workout the wavevector in the toplayer of the multilayer, if it exists.
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
		memset(subtotal,0,sizeof(subtotal));
		subtotal[0][0]=MyComplex(1,0);subtotal[1][1]=MyComplex(1,0);
		pj_mul[0] = compsqrt(qq2-*SLDmatrixREP);
	}
	
	//now calculate reflectivities
	for(ii = 0 ; ii < nlayers+1 ; ii++){
		//work out the fresnel coefficient
		if(Vmullayers>0 && ii==Vmulappend && Vmulrep>0 ){
			rj=fres(pj[ii],pj_mul[0],coefP[offset+3]);
		} else {
			rj = (ii == nlayers) ?
			((pj[ii]-pj[ii+1])/(pj[ii]+pj[ii+1]))*compexp(pj[ii]*pj[ii+1]*MyComplex(-2*subrough*subrough,0))
			:
			((pj[ii]-pj[ii+1])/(pj[ii]+pj[ii+1]))*compexp(pj[ii]*pj[ii+1]*MyComplex(-2*coefP[4*(ii+1)+7]*coefP[4*(ii+1)+7],0));
		}
		
		//work out the beta for the (non-multi)layer
		beta = (ii==0)? oneC : compexp(pj[ii] * MyComplex(0,fabs(coefP[4*ii+4])));
		
		//this is the characteristic matrix of a layer
		MI[0][0]=beta;
		MI[0][1]=rj*beta;
		MI[1][1]=oneC/beta;
		MI[1][0]=rj*MI[1][1];
		
		temp2[0][0] = MRtotal[0][0];
		temp2[0][1] = MRtotal[0][1];
		temp2[1][0] = MRtotal[1][0];
		temp2[1][1] = MRtotal[1][1];
		//multiply MR,MI to get the updated total matrix.			
		matmul(temp2,MI,MRtotal);
		
		if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0){
			//workout the wavevectors in each of the layers
			for(jj=1 ; jj < Vmullayers; jj++){
				pj_mul[jj] = compsqrt(qq2-*(SLDmatrixREP+jj));
			}
			
			//work out the fresnel coefficients
			for(jj = 0 ; jj < Vmullayers; jj++){
				rj = (jj == Vmullayers-1) ?
				//if you're in the last layer then the roughness is the roughness of the top
				((pj_mul[jj]-pj_mul[0])/(pj_mul[jj]+pj_mul[0]))*compexp((pj_mul[jj]*pj_mul[0])*MyComplex(-2*coefP[offset+3]*coefP[offset+3],0))
				:
				//otherwise it's the roughness of the layer below
				((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*MyComplex(-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3],0));
				
				
				//Beta's
				beta = compexp(MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]);
				
				MI[0][0]=beta;
				MI[0][1]=rj*beta;
				MI[1][1]=oneC/MI[0][0];
				MI[1][0]=rj*MI[1][1];
				
				temp2[0][0] = subtotal[0][0];
				temp2[0][1] = subtotal[0][1];
				temp2[1][0] = subtotal[1][0];
				temp2[1][1] = subtotal[1][1];
				
				matmul(temp2,MI,subtotal);
			};
			
			for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
				if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
					for(jj=0;jj<Vmullayers;jj++){
						beta = compexp((MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]));
						
						if(jj==Vmullayers-1){
							if(Vmulappend==nlayers){
								rj = ((pj_mul[Vmullayers-1]-pj[nlayers+1])/(pj_mul[Vmullayers-1]+pj[nlayers+1]))*compexp((pj_mul[Vmullayers-1]*pj[nlayers+1])*MyComplex(-2*subrough*subrough,0));
							} else {
								rj = ((pj_mul[Vmullayers-1]-pj[Vmulappend+1])/(pj_mul[Vmullayers-1]+pj[Vmulappend+1]))*compexp((pj_mul[Vmullayers-1]*pj[Vmulappend+1])*MyComplex(-2*coefP[4*(Vmulappend+1)+7]*coefP[4*(Vmulappend+1)+7],0));
							};
						} else {
							rj = ((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*MyComplex(-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3],0));
						}
						
						MI[0][0]=beta;
						MI[0][1]=rj*beta;
						MI[1][1]=MyComplex(1,0)/MI[0][0];
						MI[1][0]=rj*MI[1][1];
						
						temp2[0][0] = MRtotal[0][0];
						temp2[0][1] = MRtotal[0][1];
						temp2[1][0] = MRtotal[1][0];
						temp2[1][1] = MRtotal[1][1];
						
						matmul(temp2,MI,MRtotal);
					}
				} else {
					temp2[0][0] = MRtotal[0][0];
					temp2[0][1] = MRtotal[0][1];
					temp2[1][0] = MRtotal[1][0];
					temp2[1][1] = MRtotal[1][1];
					
					matmul(temp2,subtotal,MRtotal);
				};
			};
		};
		
	}
	
	den=compnorm(MRtotal[0][0]);
	
	num=compnorm(MRtotal[1][0]);
	answer=(num/den);//(num*num)/(den*den);
	answer=(answer*scale)+fabs(bkg);
	
	*result = answer;
	
done:
	if(pj != NULL)
		delete [] pj;
	if(pj_mul !=NULL)
		delete[] pj_mul;
	if(SLDmatrix != NULL)
		delete[] SLDmatrix;
	if(SLDmatrixREP != NULL)
		delete[] SLDmatrixREP;
	
	return err;	
}

