/*	Abeles.c
 A simplified project designed to act as a template for your curve fitting function.
 The fitting function is a simple polynomial. It works but is of no practical use.
 */

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#include "Abeles.h"
#include <math.h>
#include <exception>
#include "RefCalculator.h"
#ifdef _MACINTOSH_
#include <pthread.h>
#endif
#ifdef _WINDOWS_
#include "pthread.h"
#include "sched.h"
#include "semaphore.h"
#endif

/*	Abeles calculates reflectivity given a model description.
 
 Warning:
 The call to WaveData() below returns a pointer to the middle
 of an unlocked Macintosh handle. In the unlikely event that your
 calculations could cause memory to move, you should copy the coefficient
 values to local variables or an array before such operations.
 */
#include "XOPStructureAlignmentTwoByte.h"	// All structures passed to Igor are two-byte aligned.

typedef struct FitParams {
	double x;				// Independent variable.
	waveHndl waveHandle;	// Coefficient wave.
	double result;
} FitParams, *FitParamsPtr;

typedef struct FitParamsAll {
	waveHndl XWaveHandle;	// X wave (input).
	waveHndl YWaveHandle;	// Y wave (output).
	waveHndl CoefHandle;	// Coefficient wave.
	DOUBLE result;			// not actually used.
}FitParamsAll, *FitParamsAllPtr;

typedef struct SmearedParamsAll {
	waveHndl dXWaveHandle;	// X wave (input).	
	waveHndl XWaveHandle;	// X wave (input).
	waveHndl YWaveHandle;	// Y wave (output).
	waveHndl CoefHandle;	// Coefficient wave.
	DOUBLE result;			// not actually used.
}SmearedParamsAll, *SmearedParamsAllPtr;


#include "XOPStructureAlignmentReset.h"



//NUMBER OF CPUS;
int NUM_CPUS=1;

int
AbelesAll(FitParamsAllPtr p){
	long ncoefs,npoints;
	double realVal,imagVal;
	long nlayers, Vmullayers=-1, Vappendlayer=0, Vmulrep=0, err=0, ii;
	double *coefP = NULL;
	double *xP = NULL;
	double *yP = NULL;
	
	extern int NUM_CPUS;
	
	pthread_t *threads;
	refCalcParm *arg;
	long pointsEachThread = 0;
	long pointsRemaining = 0;
	long pointsConsumed = 0;
	
	if (p->CoefHandle == NULL || p->YWaveHandle == NULL || p->XWaveHandle == NULL ) 
	{
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	if (!(WaveType(p->CoefHandle) != NT_FP64 || WaveType(p->YWaveHandle) != NT_FP64 || WaveType(p->XWaveHandle) != NT_FP64
		  || WaveType(p->XWaveHandle) != NT_FP32 || WaveType(p->YWaveHandle) != NT_FP32 || WaveType(p->CoefHandle) != NT_FP32)){
		SetNaN64(&p->result);
		err = REQUIRES_SP_OR_DP_WAVE;
		goto done;
	}
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1){
		Vmullayers=0;
	} else{
		Vmullayers=(long)realVal;
	}
	
	if(FetchNumVar("Vappendlayer", &realVal, &imagVal) == -1){
		Vappendlayer = 0;
	} else{
		Vappendlayer=(long)realVal;
	}
	
	if(FetchNumVar("Vmulrep", &realVal, &imagVal) == -1){
		Vmulrep=0;
	} else{
		Vmulrep=(long)realVal;
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
	
	nlayers = (long)coefP[0];
	if(ncoefs != (4*Vmullayers+(4*nlayers+6))){
		err = INCORRECT_INPUT;
		goto done;
	};
	
	threads = (pthread_t *) malloc(NUM_CPUS * sizeof(pthread_t));
	arg=(refCalcParm *)malloc(sizeof(refCalcParm)*NUM_CPUS);
	pointsEachThread = floorl(npoints/NUM_CPUS);
	pointsRemaining = npoints;
	
	for (ii = 0; ii < NUM_CPUS; ii++){
		arg[ii].coefP = coefP;
		if(ii == NUM_CPUS-1)
			arg[ii].npoints = pointsRemaining;
		else
			arg[ii].npoints = pointsEachThread;
		
		arg[ii].Vmullayers = Vmullayers;
		arg[ii].Vappendlayer = Vappendlayer;
		arg[ii].Vmulrep = Vmulrep;
		arg[ii].xP = xP+pointsConsumed;
		arg[ii].yP = yP+pointsConsumed;
		pthread_create(&threads[ii], NULL, AbelesThreadWorker, (void *)(arg+ii));
		pointsRemaining -= pointsEachThread;
		pointsConsumed += pointsEachThread;
	}
	
	for (ii = 0; ii < NUM_CPUS; ii++){
		pthread_join(threads[ii], NULL);
	}
	
	if(threads)
		free(threads);
	if(arg)
		free(arg);
	
	if(err = MDStoreDPDataInNumericWave(p->YWaveHandle,yP))
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
	
	return err;	
}

int
smearedAbelesAll(SmearedParamsAllPtr p){
	long ncoefs,npoints;
	double realVal,imagVal;
	int nlayers,Vmullayers=-1, err=0;
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
	if (!(WaveType(p->CoefHandle) != NT_FP64 || WaveType(p->YWaveHandle) != NT_FP64 || WaveType(p->XWaveHandle) != NT_FP64
		  || WaveType(p->XWaveHandle) != NT_FP32 || WaveType(p->YWaveHandle) != NT_FP32 || WaveType(p->CoefHandle) != NT_FP32
		  || WaveType(p->dXWaveHandle) != NT_FP32 || WaveType(p->dXWaveHandle) != NT_FP64)){
		SetNaN64(&p->result);
		err = REQUIRES_SP_OR_DP_WAVE;
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
	
	if(err = smearedAbelescalcAll(coefP,yP,xP,dxP,npoints))
		goto done;
	if(err = MDStoreDPDataInNumericWave(p->YWaveHandle,yP))
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

int
Abeles_imagAll(FitParamsAllPtr p)
{
	long ncoefs,npoints;
	double *coefP = NULL;
	double *xP = NULL;
	double *yP = NULL;
	
	long nlayers, Vmullayers=-1, Vappendlayer=0, Vmulrep=0, err=0, ii;
	double realVal,imagVal;
	
	extern int NUM_CPUS;
	
	pthread_t *threads;
	refCalcParm *arg;
	long pointsEachThread = 0;
	long pointsRemaining = 0;
	long pointsConsumed = 0;

	
	if (p->CoefHandle == NIL ||	p->YWaveHandle == NIL || p->XWaveHandle == NIL ){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	if (!(WaveType(p->CoefHandle) != NT_FP64 || WaveType(p->YWaveHandle) != NT_FP64 || WaveType(p->XWaveHandle) != NT_FP64
		  || WaveType(p->XWaveHandle) != NT_FP32 || WaveType(p->YWaveHandle) != NT_FP32 || WaveType(p->CoefHandle) != NT_FP32)){
		SetNaN64(&p->result);
		err = REQUIRES_SP_OR_DP_WAVE;
		goto done;
	}
	
	ncoefs = WavePoints(p->CoefHandle);
	npoints = WavePoints(p->YWaveHandle);
	if (npoints != WavePoints(p->XWaveHandle)){
		SetNaN64(&p->result);
		err = WAVES_NOT_SAME_LENGTH;
		goto done;
	}
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1){
		Vmullayers=0;
	} else{
		Vmullayers=(long)realVal;
	}
	
	if(FetchNumVar("Vappendlayer", &realVal, &imagVal) == -1){
		Vappendlayer = 0;
	} else{
		Vappendlayer=(long)realVal;
	}
	
	if(FetchNumVar("Vmulrep", &realVal, &imagVal) == -1){
		Vmulrep=0;
	} else{
		Vmulrep=(long)realVal;
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
	
	nlayers = (long)coefP[0];
	if(ncoefs != (long)(4*Vmullayers+4*nlayers+8)){
		err = INCORRECT_INPUT;
		goto done;
	}
	
	threads = (pthread_t *) malloc(NUM_CPUS * sizeof(pthread_t));
	arg=(refCalcParm *)malloc(sizeof(refCalcParm)*NUM_CPUS);
	pointsEachThread = floorl(npoints/NUM_CPUS);
	pointsRemaining = npoints;
	
	for (ii = 0; ii < NUM_CPUS; ii++){
		arg[ii].coefP = coefP;
		if(ii == NUM_CPUS-1)
			arg[ii].npoints = pointsRemaining;
		else
			arg[ii].npoints = pointsEachThread;
		
		arg[ii].Vmullayers = Vmullayers;
		arg[ii].Vappendlayer = Vappendlayer;
		arg[ii].Vmulrep = Vmulrep;
		arg[ii].xP = xP+pointsConsumed;
		arg[ii].yP = yP+pointsConsumed;
		pthread_create(&threads[ii], NULL, AbelesImagThreadWorker, (void *)(arg+ii));
		pointsRemaining -= pointsEachThread;
		pointsConsumed += pointsEachThread;
	}
	
	for (ii = 0; ii < NUM_CPUS; ii++){
		pthread_join(threads[ii], NULL);
	}
	
	if(threads)
		free(threads);
	if(arg)
		free(arg);
	
	
	if(err = MDStoreDPDataInNumericWave(p->YWaveHandle,yP))
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
	
	return err;		
}

int
Abeles(FitParamsPtr p){
	int np,err = 0;
	double *Abelesparams = NULL;
	double x;
	char varName[MAX_OBJ_NAME+1];
	double realVal,imagVal;
	int Vmullayers;
	
	if (p->waveHandle == NULL){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	
	np= WavePoints(p->waveHandle);
	
	Abelesparams= (double*)malloc((np) * sizeof(double));	//pointer to my copy of the data.
	if(Abelesparams == NULL){
		err = NOMEM;
		goto done;
	}
	
	strcpy(varName, "Vmullayers");
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1){
		Vmullayers=0;
	} else{
		Vmullayers=(long)realVal;
	}
	
	x= p->x;
	
	if(err = MDGetDPDataFromNumericWave(p->waveHandle, Abelesparams)){
		goto done;
	}
	if((int)np!=(int)(4*Vmullayers+4*Abelesparams[0]+6)){
		err = INCORRECT_INPUT;
		goto done;
	};
	
	if(err = Abelescalc(Abelesparams,x,&p->result))
		goto done;
done:
	if(Abelesparams!=NULL)
		free(Abelesparams);
	
	return err;
}

int
Abeles_imag(FitParamsPtr p){
	int np, err = 0;
	double *Abelesparams = NULL;
	double x,result;
	char varName[MAX_OBJ_NAME+1];
	double realVal,imagVal;
	int Vmullayers;
	
	if (p->waveHandle == NULL){
		SetNaN64(&p->result);
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	
	np= WavePoints(p->waveHandle);
	
	Abelesparams= (double*)malloc((np) * sizeof(double));	//pointer to my copy of the data.
	if(Abelesparams == NULL){
		err = NOMEM;
		goto done;
	}
	
	strcpy(varName, "Vmullayers");
	if(FetchNumVar("Vmullayers", &realVal, &imagVal) == -1){
		Vmullayers=0;
	} else{
		Vmullayers=(long)realVal;
	}
	
	x= p->x;
	
	if(err = MDGetDPDataFromNumericWave(p->waveHandle, Abelesparams)){
		goto done;
	}
	if((int)np!=(int)(4*Vmullayers+4*Abelesparams[0]+8)){
		err = INCORRECT_INPUT;
		goto done;
	};
	
	if(err = Abelescalc_imag(Abelesparams,x,&result))
		goto done;
	p->result = result;
	
done:
	if(Abelesparams!=NULL)
		free(Abelesparams);
	
	return err;
}


static long
RegisterFunction()
{
	int funcIndex;
	
	funcIndex = GetXOPItem(0);			// Which function invoked ?
	switch (funcIndex) {
		case 0:							// y = Abeles(w,x) (curve fitting function).
			return((long)Abeles);	// This function is called using the direct method.
			break;
		case 1:
			return((long)Abeles_imag);
			break;
		case 2:
			return((long)AbelesAll);
			break;
		case 3:
			return((long)Abeles_imagAll);
			break;
		case 4:
			return((long)smearedAbelesAll);
			break;
	}
	return NIL;
}

/*	XOPEntry()
 
 This is the entry point from the host application to the XOP for all
 messages after the INIT message.
 */
static void
XOPEntry(void)
{	
	long result = 0;
	
	switch (GetXOPMessage()) {
		case FUNCADDRS:
			result = RegisterFunction();	// This tells Igor the address of our function.
			break;
		case CLEANUP:
#ifdef _WINDOWS_
			pthread_win32_process_detach_np();
#endif
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
#ifdef _MACINTOSH_
HOST_IMPORT int main(IORecHandle ioRecHandle)
#endif	
#ifdef _WINDOWS
HOST_IMPORT void main(IORecHandle ioRecHandle)
#endif
{	
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	extern int NUM_CPUS;

#ifdef _WINDOWS_
    SYSTEM_INFO sysInfo;  
     GetSystemInfo(&sysInfo);  
    NUM_CPUS = sysInfo.dwNumberOfProcessors;  
#endif
#ifdef _MACINTOSH_
	int mib[2];
	size_t len;
	
	mib[0] = CTL_HW;
	mib[1] = HW_NCPU;
	len = sizeof(NUM_CPUS);
	sysctl(mib, 2, &NUM_CPUS, &len, NULL, 0);
#endif

#ifdef _WINDOWS_
		pthread_win32_process_attach_np();
#endif
	
	if (igorVersion < 400)
		SetXOPResult(REQUIRES_IGOR_400);
	else
		SetXOPResult(0L);
}


MyComplex
fres(MyComplex a,MyComplex b,double rough){
	return (compexp(-2*rough*rough*a*b))*(a-b)/(a+b);
}

int
Abelescalc(double *coefP, double x, double *result){
	int err = 0;
	
	int Vmulrep=0,Vmulappend=0,Vmullayers=0;
	double realVal,imagVal;
	register int ii=0,jj=0,kk=0;
	
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
		*(SLDmatrix+ii) = 4*PI*(numtemp  - (coefP[2]*1e-6));
	}
	*(SLDmatrix) = 0;
	*(SLDmatrix+nlayers+1) = 4*PI*((coefP[3]*1e-6) - (coefP[2]*1e-6));
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal)!=-1){ // Fetch value
		Vmullayers=(int)realVal;
		if(FetchNumVar("Vappendlayer", &realVal, &imagVal)!=-1)
			Vmulappend=(int)realVal;
		if(FetchNumVar("Vmulrep", &realVal, &imagVal) !=-1) // Fetch value
			Vmulrep=(int)realVal;
		
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
			for(ii=0; ii<Vmullayers;ii+=1){
				numtemp = (coefP[3]*1e-6*coefP[(4*ii)+offset+2]/100) +(1e-6 * ((100 - coefP[(4*ii)+offset+2])/100) * coefP[(4*ii)+offset+1]);		//sld of the layer
				*(SLDmatrixREP+ii) = 4*PI*(numtemp  - (coefP[2]*1e-6));
			}
		}
	}
	
	//intialise the matrices
	memset(MRtotal,0,sizeof(MRtotal));
	MRtotal[0][0] = oneC ; MRtotal[1][1] = oneC;
	
	qq = x*x/4;
	qq2=MyComplex(qq,0);
	
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
Abelescalc_imag(double *coefP, double x, double *result){
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
		*(SLDmatrix+ii) = MyComplex(4*PI,0)*(MyComplex(coefP[4*ii+5]*1e-6,coefP[4*ii+6])-super);
	}
	*(SLDmatrix) = MyComplex(0,0);
	*(SLDmatrix+nlayers+1) = MyComplex(4*PI,0)*(sub-super);
	
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
				*(SLDmatrixREP+ii) = MyComplex(4*PI,0)*(MyComplex(coefP[(4*ii)+offset+1]*1e-6,coefP[(4*ii)+offset+2])  - super);
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


int 
smearedAbelescalcAll(double *coefP, double *yP, double *xP, double *dxP, long npoints){
	int err = 0;
	int j=0, respoints=13;
	int Vmulrep=0,Vmulappend=0,Vmullayers=0;
	double realVal=0,imagVal=0;
	int ii=0,jj=0,kk=0;
	
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
	double *dyP = NULL;
	double *ddxP = NULL;
	double *SLDmatrix = NULL;
	double *SLDmatrixREP = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		pj = new MyComplex [nlayers+2];
		SLDmatrix = new double [nlayers+2];
		dyP = new double [npoints*respoints];
		ddxP = new double [npoints*respoints];
	} catch(...){
		err = NOMEM;
		goto done;
	}
	
	memset(pj, 0, sizeof(pj));
	memset(SLDmatrix, 0, sizeof(SLDmatrix));
	memset(dyP, 0, sizeof(dyP));
	memset(ddxP, 0, sizeof(ddxP));
	
	for(ii=0 ; ii < npoints*respoints ; ii+=1){
		realVal = *(xP+ii/respoints) + (double)((ii%respoints)-(respoints-1)/2)*0.2*(*(dxP+ii/respoints));
		*(ddxP+ii) = *(xP+ii/respoints) + (double)((ii%respoints)-(respoints-1)/2)*0.2*(*(dxP+ii/respoints));
	}
	
	scale = coefP[1];
	bkg = fabs(coefP[4]);
	subrough = coefP[5];
	
	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 6;
	
	//fillout all the SLD's for all the layers
	for(ii=1; ii<nlayers+1;ii+=1){
		numtemp = 1.e-6 * ((100. - coefP[4*ii+4])/100.) * coefP[4*ii+3]+ (coefP[4*ii+4]*coefP[3]*1.e-6)/100.;		//sld of the layer
		*(SLDmatrix+ii) = 4*PI*(numtemp  - (coefP[2]*1e-6));
	}
	*(SLDmatrix) = 0;
	*(SLDmatrix+nlayers+1) = 4*PI*((coefP[3]*1e-6) - (coefP[2]*1e-6));
	
	if(FetchNumVar("Vmullayers", &realVal, &imagVal)!=-1){ // Fetch value
		Vmullayers=(int)realVal;
		if(FetchNumVar("Vappendlayer", &realVal, &imagVal)!=-1) // Fetch value
			Vmulappend=(int)realVal;
		if(FetchNumVar("Vmulrep", &realVal, &imagVal) !=-1) // Fetch value
			Vmulrep=(int)realVal;
		
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
			for(ii=0; ii<Vmullayers;ii+=1){
				numtemp = (coefP[3]*1e-6*coefP[(4*ii)+offset+2]/100) +(1e-6 * ((100 - coefP[(4*ii)+offset+2])/100) * coefP[(4*ii)+offset+1]);		//sld of the layer
				*(SLDmatrixREP+ii) = 4*PI*(numtemp  - (coefP[2]*1e-6));
			}
		}
	}
	
	for (j = 0; j < npoints * respoints; j++) {
		//intialise the matrices
		memset(MRtotal,0,sizeof(MRtotal));
		MRtotal[0][0] = oneC ; MRtotal[1][1] = oneC;
		
		qq = ddxP[j]*ddxP[j]/4;
		qq2=MyComplex(qq,0);
		
		for(ii=0; ii<nlayers+2 ; ii++){			//work out the wavevector in each of the layers
			pj[ii] = (*(SLDmatrix+ii)>qq) ? compsqrt(qq2-MyComplex(*(SLDmatrix+ii),0)): MyComplex(sqrt(qq-*(SLDmatrix+ii)),0);
		}
		
		//workout the wavevector in the toplayer of the multilayer, if it exists.
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
			memset(subtotal,0,sizeof(subtotal));
			subtotal[0][0]=MyComplex(1,0);subtotal[1][1]=MyComplex(1,0);
			pj_mul[0] = (*(SLDmatrixREP)>qq) ? compsqrt(qq2-MyComplex(*SLDmatrixREP,0)): MyComplex(sqrt(qq-*SLDmatrixREP),0);
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
		answer=((num/den)*scale)+bkg;
		
		*(dyP+j) = answer;
	}
	
	for(ii=0 ; ii < npoints ; ii+=1){
		*(yP+ii) = 0.056*(*(dyP+ii*respoints));
		*(yP+ii) += 0.135*(*(dyP+ii*respoints+1));														
		*(yP+ii) += 0.278*(*(dyP+ii*respoints+2));
		*(yP+ii) += 0.487*(*(dyP+ii*respoints+3));
		*(yP+ii) += 0.726*(*(dyP+ii*respoints+4));
		*(yP+ii) += 0.923*(*(dyP+ii*respoints+5));
		*(yP+ii) += (*(dyP+ii*respoints+6));
		*(yP+ii) += 0.923*(*(dyP+ii*respoints+7));
		*(yP+ii) += 0.726*(*(dyP+ii*respoints+8));
		*(yP+ii) += 0.487*(*(dyP+ii*respoints+9));
		*(yP+ii) += 0.278*(*(dyP+ii*respoints+10));
		*(yP+ii) += 0.135*(*(dyP+ii*respoints+11));
		*(yP+ii) += 0.056*(*(dyP+ii*respoints+12));
		*(yP+ii) /=6.211;
	}
	
done:
	if(pj != NULL)
		delete [] pj;
	if(pj_mul !=NULL)
		delete[] pj_mul;
	if(SLDmatrix != NULL)
		delete[] SLDmatrix;
	if(SLDmatrixREP != NULL)
		delete[] SLDmatrixREP;
	if(ddxP != NULL)
		delete[] ddxP;
	if(dyP != NULL)
		delete[] dyP;
	
	return err;
}