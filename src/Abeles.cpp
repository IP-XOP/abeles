/*	Abeles.cpp
calculates specular reflectivity as a function of Q momentum transfer.
 
 Copyright Andrew Nelson and ANSTO 2007
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#include "Abeles.h"
#include <math.h>
#include <exception>
#include "RefCalculator.h"
#include <vector>


void matmul(complex<double> a[2][2], complex<double> b[2][2], complex<double> c[2][2]){
    c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0];
    c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1];
    c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0];
    c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1];
}

/*
 because we are going to do the calculation in a threaded fashion we need to know the number of CPU's.
 This is stored in a global so we don't need to find it out every time we call the function.
 */
int NUM_CPUS = 1;
GaussWeights *pinstance = NULL;
pthread_mutex_t changeWeightsMutex = PTHREAD_MUTEX_INITIALIZER;

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

#pragma pack()	// All structures passed to Igor are two-byte aligned.



void getMultiLayerParams(int *Vmullayers, int *Vappendlayer, int *Vmulrep){
	//temporary places for calculation
	double realVal, imagVal;
	
	//this may meanthere is a multilayer model specified
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1)
		*Vmullayers = 0;
	else
		*Vmullayers=(int)realVal;
	
	if(FetchNumVar("Vappendlayer", &realVal, &imagVal) == -1)
		*Vappendlayer = 0;
	else
		*Vappendlayer=(int)realVal;
	
	if(FetchNumVar("Vmulrep", &realVal, &imagVal) == -1)
		*Vmulrep=0;
	else
		*Vmulrep=(int)realVal;
	
}


/*mode:
0 = AbelesAll
1 = Abeles_imagAll
 */
 
int AbelesAllWrapper(FitParamsAllPtr p, int mode){
	int err = 0;

	//number of coefficients, number of Q points
	CountInt ncoefs,npoints, smearedPoints;
	
	//how many layers you have
	long nlayers;
	
	//how many layers in multilayer (if you have one)
	int Vmullayers =- 1;
	//where in the normal model the multilayer gets appended to
	int Vappendlayer=0;
	//how many times the multilayer repeats
	int Vmulrep=0;
	
	//a variable for iterating for loops
	CountInt ii;
	
	//we will extract values from the supplied waves and store them as double arrays
	//these pointers refer to the values extracted.
	double *coefP = NULL;
    double bkd = 0.0;
	double *xP = NULL;
	double *yP = NULL;
	double *dxP = NULL;
    double *dp = NULL;
	double *calcX, *calcY;
	double *smearedX = NULL;
	double *smearedY = NULL;
    double real = 0;
    double imag = 0;
    int retval = -1;
	
	//variables for the threadwise calculation of the reflectivity
	int threadsToCreate = 1;
	double isItWorthThreading = 0;
	pthread_t *threads = NULL;
	refCalcParm *arg = NULL;
	CountInt pointsEachThread = 0;
	CountInt pointsRemaining = 0;
	CountInt pointsConsumed = 0;
	
	//values for the gaussian quadrature.
	int isSmeared = 0;
	int RESPOINTS = 17;
    const double m_pi = 3.14159265358979323846;
	const double INTLIMIT = 3.5;		//integration between -3.5 and 3.5 sigma
	const double FWHM = 2 * sqrt(2 * log(2.0));
	double va, vb, sigma;
    std::vector<double> weights;
    std::vector<double> abscissa;
	
	//the functions that will do the calculation
	void* (*threadWorkerFunc)(void*);
	int (*calcAllFunc)(const double *, double *, const double *,long, int, int, int);
			
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
	if(npoints == WavePoints(p->XWaveHandle)){
		isSmeared = 0;
		smearedPoints = npoints;
	}
	else if(npoints * 2 == WavePoints(p->XWaveHandle)){
		isSmeared = 1;
        retval = FetchNumVar("V_gausspoints", &real, &imag);
        if(retval != -1)
            RESPOINTS = (int) real;
		smearedPoints = npoints * RESPOINTS;        
	} else {
		err = WAVES_NOT_SAME_LENGTH;
	}
	
	coefP = (double*) WaveData(p->CoefHandle);	
	// the number of layers should relate to the size of the coefficient wave (4*N+6 + 4*Vmullayers)
	// if not, the coefficient wave is wrong.
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
            bkd = coefP[4];			
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
            bkd = coefP[6];            
			break;
		default:
			break;
	}
	
	if(isSmeared){
		smearedX = (double*) malloc(sizeof(double) * RESPOINTS * npoints);
		if(!smearedX){
			err = NOMEM;
			goto done;
		}
		smearedY = (double*) malloc(sizeof(double) * smearedPoints);
		if(!smearedY){
			err = NOMEM;
			goto done;
		}
		yP = (double*) WaveData(p->YWaveHandle);
		xP = (double*) WaveData(p->XWaveHandle);
		dxP = xP + npoints;
		
        pinstance->getGaussWeight(RESPOINTS, abscissa, weights);
        
		for(ii = 0 ; ii < smearedPoints ; ii += 1){
            CountInt idx = ii / RESPOINTS;
            sigma = dxP[idx] / FWHM;
            va = -INTLIMIT * sigma + xP[idx];
            vb = INTLIMIT * sigma + xP[idx];
			smearedX[ii] = (abscissa[ii % RESPOINTS] * (vb-va) + vb + va)/2;
		}
        
        for(ii = 0 ; ii < RESPOINTS ; ii++)
            weights[ii] = exp(-0.5 * pow(INTLIMIT * abscissa[ii], 2)) * weights[ii] / sqrt(2 * m_pi);

		calcY = smearedY;
		calcX = smearedX;
        switch(mode){
            case 0:
                coefP[4] = 0;
                break;
            case 1:
                coefP[6] = 0;
                break;
        }
	} else {
		yP = (double*) WaveData(p->YWaveHandle);
		xP = (double*) WaveData(p->XWaveHandle);
		
		calcX = xP;
		calcY = yP;
	}
		
	//this relationship was worked out for a dualcore machine.
	//I worked out how long it took for a certain number of points and a certain number of layers
	//then I calculated the cross over for when it was worth threading, rather than not.
	//i.e. for a given number of layers, how many points were required for multithreading to be worthwhile.
	//I plotted the number of points (y) vs the number of layers (x), giving the following relationship.
	isItWorthThreading = 3.382 + 641. *pow(coefP[0], -0.73547);
	if((double) smearedPoints < isItWorthThreading)
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
	pointsEachThread = (CountInt) floorl((long) smearedPoints / threadsToCreate);
	pointsRemaining = smearedPoints;
	
	for (ii = 0; ii < threadsToCreate - 1; ii++){
		arg[ii].coefP = (double*) WaveData(p->CoefHandle);
		arg[ii].npoints = (long) pointsEachThread;
		
		arg[ii].Vmullayers = Vmullayers;
		arg[ii].Vappendlayer = Vappendlayer;
		arg[ii].Vmulrep = Vmulrep;
		
		//the following two lines specify where the Q values and R values will be sourced/written.
		//i.e. an offset of the original array.
		arg[ii].xP = calcX + pointsConsumed;
		arg[ii].yP = calcY + pointsConsumed;
		
		pthread_create(&threads[ii], NULL, threadWorkerFunc, (void *)(arg+ii));
		pointsRemaining -= pointsEachThread;
		pointsConsumed += pointsEachThread;
	}
	
	//do the remaining points in the main thread.
	if(err = calcAllFunc((double*) WaveData(p->CoefHandle), calcY+pointsConsumed, calcX+pointsConsumed, (long) pointsRemaining, Vmullayers, Vappendlayer, Vmulrep))
		goto done;
	
	for (ii = 0; ii < threadsToCreate - 1 ; ii++)
		pthread_join(threads[ii], NULL);
	
	/**if you're not smeared, you are all finished.  If you aren't, then you still have to calculate the smearing.
	and refill the wave
	 **/
	if(isSmeared){
        memset(yP, 0, sizeof(double) * npoints);
        switch(mode){
            case 0:
                coefP[4] = bkd;
                break;
            case 1:
                coefP[6] = bkd;
                break;
        }
        
		for(ii = 0 ; ii < smearedPoints ; ii += 1)
            yP[ii / RESPOINTS] += calcY[ii] * weights[ii % RESPOINTS];
            
        dp = yP;
        for(ii = 0 ; ii < npoints ; ii++, dp++)
            *dp = (*dp * INTLIMIT) + bkd;
	}
	
	WaveHandleModified(p->YWaveHandle);
	p->result = 0;		// not actually used by FuncFit
	
done:
	if(threads)
		free(threads);
	if(arg)
		free(arg);
	if(smearedX)
		free(smearedX);
	if(smearedY)
		free(smearedY);
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
smearedAbelesAll(FitParamsAllPtr p){
	return AbelesAllWrapper(p, 0);	
}


extern "C" int
Abeles(FitParamsPtr p){
	int err = 0;
	CountInt np = 0;
	double *Abelesparams = NULL;
	double x;
	int Vmullayers = 0;
	
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

	Vmullayers = (((int)np - 6) / 4) - (int) Abelesparams[0];
	
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
	int err = 0;
	CountInt np = 0;
	double *Abelesparams = NULL;
	double x,result;
	int Vmullayers = 0;
	
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
	
	Vmullayers = (((int)np - 8) / 4) - (int) Abelesparams[0];

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


static XOPIORecResult
RegisterFunction()
{
	XOPIORecResult funcIndex;
	
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
        case CLEANUP:
            delete(pinstance);
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
HOST_IMPORT int XOPMain(IORecHandle ioRecHandle){
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
    
    pinstance = new GaussWeights();
    if(!pinstance)
        return NOMEM;

	// find out the number of CPU's.
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
	
	if (igorVersion < 700){
		SetXOPResult(IGOR_OBSOLETE);
		return EXIT_FAILURE;
	}else{
		SetXOPResult(0L);
		return EXIT_SUCCESS;
	}
}


complex<double> fres(complex<double> a, complex<double> b, double rough){
    return (std::exp(-2. * rough * rough * a * b)) * (a - b) / (a + b);
}

int
Abelescalc(const double *coefP, double x, double *result){
	int err = 0;
	
	int Vmulrep = 0,Vmulappend = 0, Vmullayers=0;
	double realVal, imagVal;
	register int ii=0, jj=0, kk=0;
	
	double scale,bkg,subrough;
	double num=0,den=0, answer=0,qq;
	complex<double> temp,SLD,beta,rj;
	double numtemp=0;
	int offset=0;
	complex<double>  MRtotal[2][2];
	complex<double> subtotal[2][2];
	complex<double> MI[2][2];
	complex<double> temp2[2][2];
	complex<double> qq2;
	complex<double> oneC = complex<double>(1, 0);
	complex<double> *pj_mul = NULL;
	complex<double> *pj = NULL;
	double *SLDmatrix = NULL;
	double *SLDmatrixREP = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		pj = new complex<double> [nlayers+2];
		SLDmatrix = new double [nlayers+2];
	} catch(...){
		err = NOMEM;
		goto done;
	}
	
	memset(pj, 0, (nlayers + 2) * sizeof(*pj));
	memset(SLDmatrix, 0, (nlayers + 2) * sizeof(*SLDmatrix));
	
	scale = coefP[1];
	bkg = (coefP[4]);
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
				pj_mul = new complex<double> [Vmullayers];
			} catch(...){
				err = NOMEM;
				goto done;
			}
			memset(pj_mul, 0, Vmullayers * sizeof(*pj_mul));
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
	qq2 = complex<double>(qq, 0);
	
	for(ii=0; ii<nlayers+2 ; ii++)		//work out the wavevector in each of the layers
        pj[ii] = std::sqrt(qq2 - *(SLDmatrix + ii));
	
	// workout the wavevector in the toplayer of the multilayer, if it exists.
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
		memset(subtotal,0,sizeof(subtotal));
		subtotal[0][0] = oneC;
        subtotal[1][1] = oneC;
        pj_mul[0] = std::sqrt(qq2 - *SLDmatrixREP);
	}
	
	//now calculate reflectivities
	for(ii = 0 ; ii < nlayers+1 ; ii++){
		//work out the fresnel coefficients
		//this looks more complicated than it really is.
		//the reason it looks so convoluted is because if there is no complex part of the wavevector,
		//then it is faster to do the calc with real arithmetic then put it into a complex number.
		if(Vmullayers>0 && ii==Vmulappend && Vmulrep>0 ){
			rj = fres(pj[ii], pj_mul[0], coefP[offset+3]);
		} else {
            rj = (ii == nlayers) ?
            fres(pj[ii], pj[ii+1], subrough)
            :
            fres(pj[ii], pj[ii+1], coefP[4*(ii+1)+5]);
		}
		
		//work out the beta for the (non-multi)layer
        beta = (ii==0)? oneC : std::exp(pj[ii] * complex<double>(0, fabs(coefP[4*ii+2])));
		
		//this is the characteristic matrix of a layer
		MI[0][0] = beta;
		MI[0][1] = rj*beta;
		MI[1][1] = oneC / beta;
		MI[1][0]=rj * MI[1][1];
		
		temp2[0][0] = MRtotal[0][0];
		temp2[0][1] = MRtotal[0][1];
		temp2[1][0] = MRtotal[1][0];
		temp2[1][1] = MRtotal[1][1];
		//multiply MR,MI to get the updated total matrix.			
		matmul(temp2, MI, MRtotal);
		
		if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0){
			//workout the wavevectors in each of the layers
			for(jj=1 ; jj < Vmullayers; jj++){
                pj_mul[jj] = std::sqrt(qq2 - *(SLDmatrixREP + jj));
			}
			
			//work out the fresnel coefficients
			for(jj = 0 ; jj < Vmullayers; jj++){
				
				rj = (jj == Vmullayers-1) ?
				//if you're in the last layer then the roughness is the roughness of the top
                fres(pj_mul[jj], pj_mul[0], coefP[offset+3])
				:
				//otherwise it's the roughness of the layer below
                fres(pj_mul[jj], pj_mul[jj+1], coefP[4*(jj+1)+offset+3]);

				// Beta's
                beta = std::exp(complex<double>(0, fabs(coefP[4*jj+offset])) * pj_mul[jj]);
				
				MI[0][0] = beta;
				MI[0][1] = rj * beta;
				MI[1][1] = oneC / beta;
				MI[1][0] = rj * MI[1][1];
				
				temp2[0][0] = subtotal[0][0];
				temp2[0][1] = subtotal[0][1];
				temp2[1][0] = subtotal[1][0];
				temp2[1][1] = subtotal[1][1];
				
				matmul(temp2,MI,subtotal);
			};
			
			for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
				if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
					for(jj=0;jj<Vmullayers;jj++){
                        beta = std::exp(complex<double>(0, fabs(coefP[4*jj+offset])) * pj_mul[jj]);
						
						if(jj==Vmullayers-1){
							if(Vmulappend==nlayers){
                                rj = ((pj_mul[Vmullayers-1]-pj[nlayers+1])/(pj_mul[Vmullayers-1]+pj[nlayers+1])) * std::exp(pj_mul[Vmullayers-1]*pj[nlayers+1] * -2. * subrough * subrough);
							} else {
                                rj = ((pj_mul[Vmullayers-1]-pj[Vmulappend+1])/(pj_mul[Vmullayers-1]+pj[Vmulappend+1]))*std::exp(pj_mul[Vmullayers-1] * pj[Vmulappend+1] * -2. * coefP[4*(Vmulappend+1)+5] * coefP[4*(Vmulappend+1)+5]);
							};
						} else {
                            rj = ((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1])) * std::exp(pj_mul[jj] * pj_mul[jj+1] * -2. * coefP[4*(jj+1)+offset+3] * coefP[4*(jj+1)+offset+3]);
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
	
    num = std::norm(MRtotal[1][0]);
    den = std::norm(MRtotal[0][0]);
    answer = (num / den);
    answer = (answer * scale) + bkg;
	
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
	
	complex<double> super;
	complex<double> sub;
	complex<double> temp,SLD,beta,rj,arg;
	complex<double> oneC = complex<double>(1,0);
	int offset=0;
	complex<double> MRtotal[2][2];
	complex<double> subtotal[2][2];
	complex<double> MI[2][2];
	complex<double> temp2[2][2];
	complex<double> qq2;
	complex<double> *pj_mul = NULL;
	complex<double> *pj = NULL;
	complex<double> *SLDmatrix = NULL;
	complex<double> *SLDmatrixREP = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		pj = new complex<double>[nlayers+2];
		SLDmatrix = new complex<double>[nlayers+2];
	} catch(...){
		err = NOMEM;
		goto done;
	}
	
	memset(pj, 0, (nlayers + 2) * sizeof(*pj));
	memset(SLDmatrix, 0, (nlayers + 2) * sizeof(*SLDmatrix));
	
	scale = coefP[1];
	bkg = coefP[6];
	subrough = coefP[7];
	sub= complex<double>(coefP[4]*1e-6, coefP[5]);
	super = complex<double>(coefP[2]*1e-6, 0);
	
	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 8;
	
	//fillout all the SLD's for all the layers
	for(ii=1; ii<nlayers+1;ii+=1)
		*(SLDmatrix+ii) = 4. * PI * (complex<double>(coefP[4*ii+5] * 1e-6, coefP[4*ii+6] + 1.e-100) - super);

	*(SLDmatrix) = complex<double>(0, 0);
	*(SLDmatrix+nlayers + 1) = 4. * PI * (sub - super);
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal)!=-1){ // Fetch value
		Vmullayers=(int)realVal;
		if(FetchNumVar("Vappendlayer", &realVal, &imagVal)!=-1) // Fetch value
			Vmulappend=(int)realVal;
		if(FetchNumVar("Vmulrep", &realVal, &imagVal) !=-1) // Fetch value
			Vmulrep=(int)realVal;
		
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
			//set up an array for wavevectors
			try{
				SLDmatrixREP = new complex<double>[Vmullayers];
				pj_mul = new complex<double>[Vmullayers];
			} catch(...){
				err = NOMEM;
				goto done;
			}
			memset(pj_mul, 0, Vmullayers * sizeof(*pj_mul));
			memset(SLDmatrixREP, 0, Vmullayers * sizeof(*SLDmatrixREP));
			for(ii=0; ii<Vmullayers;ii+=1){
				*(SLDmatrixREP+ii) = 4 * PI * (complex<double>(coefP[(4*ii)+offset+1] * 1e-6, coefP[(4*ii)+offset+2] + 1.e-100)  - super);
			}
		}
	}
	
	//intialise the matrices
	memset(MRtotal, 0, sizeof(MRtotal));
	MRtotal[0][0] = oneC;
    MRtotal[1][1] = oneC;
	
	qq2 = complex<double>(x * x / 4, 0);
	
	for(ii=0; ii<nlayers+2 ; ii++){			//work out the wavevector in each of the layers
        pj[ii] = std::sqrt(qq2 - *(SLDmatrix+ii));
	}
	
	//workout the wavevector in the toplayer of the multilayer, if it exists.
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
		memset(subtotal,0,sizeof(subtotal));
		subtotal[0][0] = oneC;
        subtotal[1][1] = oneC;
        pj_mul[0] = std::sqrt(qq2 - *SLDmatrixREP);
	}
	
	//now calculate reflectivities
	for(ii = 0 ; ii < nlayers+1 ; ii++){
		//work out the fresnel coefficient
		if(Vmullayers>0 && ii==Vmulappend && Vmulrep>0 ){
			rj=fres(pj[ii], pj_mul[0], coefP[offset+3]);
		} else {
			rj = (ii == nlayers) ?
			fres(pj[ii], pj[ii+1], subrough)
			:
            fres(pj[ii], pj[ii+1], coefP[4*(ii+1) + 7]);
        }
		
		//work out the beta for the (non-multi)layer
        beta = (ii==0)? oneC : std::exp(pj[ii] * complex<double>(0, fabs(coefP[4*ii+4])));
		
		//this is the characteristic matrix of a layer
		MI[0][0] = beta;
		MI[0][1] = rj * beta;
		MI[1][1] = oneC / beta;
		MI[1][0] = rj * MI[1][1];
		
		temp2[0][0] = MRtotal[0][0];
		temp2[0][1] = MRtotal[0][1];
		temp2[1][0] = MRtotal[1][0];
		temp2[1][1] = MRtotal[1][1];
		//multiply MR,MI to get the updated total matrix.			
		matmul(temp2, MI, MRtotal);
		
		if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0){
			//workout the wavevectors in each of the layers
			for(jj=1 ; jj < Vmullayers; jj++){
                pj_mul[jj] = std::sqrt(qq2 - *(SLDmatrixREP + jj));
			}
			
			//work out the fresnel coefficients
			for(jj = 0 ; jj < Vmullayers; jj++){
				rj = (jj == Vmullayers-1) ?
				//if you're in the last layer then the roughness is the roughness of the top
                fres(pj_mul[jj], pj_mul[0], coefP[offset+3])
				:
                //otherwise it's the roughness of the layer below
                fres(pj_mul[jj], pj_mul[jj+1], coefP[4*(jj+1)+offset+3]);

				//Beta's
                beta = std::exp(complex<double>(0, fabs(coefP[4*jj+offset])) * pj_mul[jj]);
				
				MI[0][0] = beta;
				MI[0][1] = rj * beta;
				MI[1][1] = oneC / MI[0][0];
				MI[1][0] = rj * MI[1][1];
				
				temp2[0][0] = subtotal[0][0];
				temp2[0][1] = subtotal[0][1];
				temp2[1][0] = subtotal[1][0];
				temp2[1][1] = subtotal[1][1];
				
				matmul(temp2, MI, subtotal);
			};
			
			for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
				if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
					for(jj=0;jj<Vmullayers;jj++){
                        beta = std::exp(complex<double>(0, fabs(coefP[4*jj+offset])) * pj_mul[jj]);
						
						if(jj==Vmullayers-1){
							if(Vmulappend==nlayers){
                                rj = fres(pj_mul[Vmullayers-1], pj[nlayers+1], subrough);
							} else {
                                rj = fres(pj_mul[Vmullayers-1], pj[Vmulappend+1], coefP[4*(Vmulappend+1)+7]);
							};
						} else {
                            rj = fres(pj_mul[jj], pj_mul[jj+1], coefP[4*(jj+1)+offset+3]);
						}
						
						MI[0][0] = beta;
						MI[0][1] = rj * beta;
						MI[1][1] = oneC / MI[0][0];
						MI[1][0] = rj*MI[1][1];
						
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
					
					matmul(temp2, subtotal, MRtotal);
				};
			};
		};
		
	}
	
    num = std::norm(MRtotal[1][0]);
    den = std::norm(MRtotal[0][0]);
    answer = (num / den);
    answer = (answer * scale) + bkg;
	
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
