/*
 *  Threadcalculator.cpp
 *  Abeles
 *
 Copyright Andrew Nelson and ANSTO 2009
 *
 */
#include "XOPStandardHeaders.h"
#include "RefCalculator.h"
#include "MyComplex.h"
#include "Abeles.h"

#ifdef MACIGOR
	#include <pthread.h>
#endif

#ifdef WINIGOR
#include "pthread.h"
#include "sched.h"
#include "semaphore.h"
#endif

using namespace MyComplexNumber;

void *AbelesThreadWorker(void *arg){
	int err = NULL;
	refCalcParm *p = (refCalcParm *) arg;
	err = AbelesCalcAll(p->coefP, p->yP, p->xP, p->npoints, p->Vmullayers, p->Vappendlayer, p->Vmulrep);
//	err = ParrattCalcAll(p->coefP, p->yP, p->xP, p->npoints);
	
	pthread_exit((void*)err);
	return NULL;
}

void *AbelesImagThreadWorker(void *arg){
	int err = NULL;
	refCalcParm *p = (refCalcParm *) arg;
	err = AbelesCalc_ImagAll(p->coefP, p->yP, p->xP, p->npoints, p->Vmullayers, p->Vappendlayer, p->Vmulrep);
	
	pthread_exit((void*)err);
	return NULL;
}

void *realReflectanceThreadWorker(void *arg){
	int err = NULL;
	refCalcParm *p = (refCalcParm *) arg;
	err = realReflectance(p->coefP, p->yP, p->xP, p->npoints);
	
	pthread_exit((void*)err);
	return NULL;
}

int ParrattCalcAll(const double *coefP, double *yP, const double *xP,long npoints){
	int err = 0;
	int j;
	
	int ii=0;
	double scale,bkg,subrough;
	double answer=0, qq;
	double anum,anum2;
	MyComplex temp, beta, rj, RRJ, RRJ_1, kzj, kzj_1;
	double numtemp=0;
	int offset=0;

	double *SLDmatrix = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		SLDmatrix = new double [nlayers+2];
	} catch(...){
		err = NOMEM;
		goto done;
	}
	scale = coefP[1];
	bkg = fabs(coefP[4]);
	subrough = coefP[5];
	
	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 6;
	
	//fillout all the SLD's for all the layers
	for(ii = 1 ; ii < nlayers+1 ; ii += 1){
		numtemp = 1.e-6 * ((100. - coefP[4*ii+4])/100.) * coefP[4*ii+3]+ (coefP[4*ii+4]*coefP[3]*1.e-6)/100.;		//sld of the layer
		*(SLDmatrix+ii) = 4 * PI * (numtemp);
	}
	*(SLDmatrix) = (coefP[2] * 1e-6);
	*(SLDmatrix+nlayers+1) = 4 * PI * ((coefP[3] * 1e-6));
	
	for (j = 0; j < npoints ; j+=1) {
		//intialise the matrices
		
		qq = xP[j]*xP[j]/4;

		//start from subphase
		kzj_1 = (*(SLDmatrix+nlayers+1) > qq) ? MyComplex(0, sqrt(fabs(qq - *(SLDmatrix+nlayers+1)))) : MyComplex(sqrt(qq - *(SLDmatrix+nlayers+1)), 0);		
		RRJ_1.re = 0;
		RRJ_1.im = 0;
		
		for(ii = nlayers ; ii > -1 ; ii--){
			//wave vector in layer ii
//			kzj = (*(SLDmatrix+ii) > qq) ? MyComplex(0, sqrt(fabs(qq - *(SLDmatrix+ii)))) : MyComplex(sqrt(qq - *(SLDmatrix+ii)), 0);
			if((*(SLDmatrix+ii) > qq)) {
				kzj.re = 0;
				kzj.im = sqrt(fabs(qq - *(SLDmatrix+ii)));
			} else {
				kzj.re = sqrt(qq - *(SLDmatrix+ii));
				kzj.im = 0;
			}

			//work out reflection coefficient, rj
			if(kzj.im == 0 && kzj_1.im == 0){
				anum = kzj.re;
				anum2 = kzj_1.re;
				rj.re = (ii==nlayers) ? 
				((anum-anum2) / (anum+anum2)) * exp(anum*anum2*-2*subrough*subrough)
				:
				((anum-anum2) / (anum+anum2)) * exp(anum*anum2*-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5]);
				rj.im = 0.;
			} else {
				rj = (ii == nlayers) ?
				((kzj - kzj_1) / (kzj + kzj_1)) * compexp(kzj * kzj_1 * -2 * subrough * subrough)
				:
				((kzj-kzj_1)/(kzj+kzj_1))*compexp(kzj*kzj_1*-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5]);	
			};
			
			if(ii == nlayers){
				RRJ.re = rj.re;
				RRJ.im = rj.im;
			} else {
				temp.re = - 2 * fabs(coefP[4 * (ii + 1) + 2]) * kzj_1.im ;
				temp.im = kzj_1.re * -2 *  fabs(coefP[4 * (ii + 1) + 2]);
				beta = compexp(temp);
				RRJ = (rj + RRJ_1 * beta) / (1 + rj * RRJ_1 * beta);
		/*		double numRE,numIM,denomRE,denomIM, denom2;
				numRE =  RRJ_1.re*beta.re - RRJ_1.im*beta.im;
				numIM =  RRJ_1.re*beta.im + RRJ_1.im*beta.re;
				
				denomRE = 1+ numRE * rj.re - numIM * rj.im;
				denomIM = numRE * rj.im + numIM * rj.re;
				
				numRE += rj.re;
				numIM += rj.im;
				
				denom2 = denomRE * denomRE + denomIM * denomIM;
				RRJ.re = (numRE * denomRE + numIM * denomIM) / denom2;
				RRJ.im = (numIM * denomRE - numRE * denomIM) / denom2;		*/
			}
			
			kzj_1.re = kzj.re;
			kzj_1.im = kzj.im;
			RRJ_1.re = RRJ.re;
			RRJ_1.im = RRJ.im;
		}
		answer = (scale*(RRJ.im * RRJ.im + RRJ.re * RRJ.re)) + bkg; 
				
		*yP++ = answer;
	}
	
done:
	if(SLDmatrix != NULL)
		delete[] SLDmatrix;
	
	return err;
}

int realReflectance(const double *coefP, double *yP, const double *xP,long npoints){
	int err = 0;
	int j;
	
	int ii=0;
	double scale,bkg,subrough;
	double answer=0, qq;
	double anum,anum2;
	MyComplex temp, beta, rj, RRJ, RRJ_1, kzj, kzj_1;
	double numtemp=0;
	int offset=0;
	
	double *SLDmatrix = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		SLDmatrix = new double [nlayers+2];
	} catch(...){
		err = NOMEM;
		goto done;
	}
	scale = coefP[1];
	bkg = fabs(coefP[4]);
	subrough = coefP[5];
	
	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 6;
	
	//fillout all the SLD's for all the layers
	for(ii = 1 ; ii < nlayers+1 ; ii += 1){
		numtemp = 1.e-6 * ((100. - coefP[4*ii+4])/100.) * coefP[4*ii+3]+ (coefP[4*ii+4]*coefP[3]*1.e-6)/100.;		//sld of the layer
		*(SLDmatrix+ii) = 4 * PI * (numtemp  - (coefP[2] * 1e-6));		
	}
	*(SLDmatrix) = 0;
	*(SLDmatrix+nlayers+1) = 4 * PI * ((coefP[3] * 1e-6) - (coefP[2] * 1e-6));
	
	for (j = 0; j < npoints ; j+=1) {
		//intialise the matrices
		
		qq = xP[j]*xP[j]/4;
		
		//start from subphase
		kzj_1 = (*(SLDmatrix+nlayers+1) > qq) ? MyComplex(0, sqrt(fabs(qq - *(SLDmatrix+nlayers+1)))) : MyComplex(sqrt(qq - *(SLDmatrix+nlayers+1)), 0);		
		RRJ_1.re = 0;
		RRJ_1.im = 0;
		
		for(ii = nlayers ; ii > -1 ; ii--){
			//wave vector in layer ii
			if((*(SLDmatrix+ii) > qq)) {
				kzj.re = 0;
				kzj.im = sqrt(fabs(qq - *(SLDmatrix+ii)));
			} else {
				kzj.re = sqrt(qq - *(SLDmatrix+ii));
				kzj.im = 0;
			}
			
			//work out reflection coefficient, rj
			if(kzj.im == 0 && kzj_1.im == 0){
				anum = kzj.re;
				anum2 = kzj_1.re;
				rj.re = (ii==nlayers) ? 
				((anum-anum2) / (anum+anum2)) * exp(anum*anum2*-2*subrough*subrough)
				:
				((anum-anum2) / (anum+anum2)) * exp(anum*anum2*-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5]);
				rj.im = 0.;
			} else {
				rj = (ii == nlayers) ?
				((kzj - kzj_1) / (kzj + kzj_1)) * compexp(kzj * kzj_1 * -2 * subrough * subrough)
				:
				((kzj-kzj_1)/(kzj+kzj_1))*compexp(kzj*kzj_1*-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5]);	
			};
			
			if(ii == nlayers){
				RRJ.re = rj.re;
				RRJ.im = rj.im;
			} else {
				temp.re = - 2 * fabs(coefP[4 * (ii + 1) + 2]) * kzj_1.im ;
				temp.im = kzj_1.re * -2 *  fabs(coefP[4 * (ii + 1) + 2]);
				beta = compexp(temp);
				RRJ = (rj + RRJ_1 * beta) / (1 + rj * RRJ_1 * beta);

			}
			
			kzj_1.re = kzj.re;
			kzj_1.im = kzj.im;
			RRJ_1.re = RRJ.re;
			RRJ_1.im = RRJ.im;
		}
		
		answer = RRJ.re;
		
		*yP++ = answer;
	}
	
done:
	if(SLDmatrix != NULL)
		delete[] SLDmatrix;
	
	return err;
}

int 
AbelesCalcAll(const double *coefP, double *yP, const double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep){
	int err = 0;
	int j;
	
	int ii=0,jj=0,kk=0;

	double scale,bkg,subrough;
	double num=0, den=0, answer=0, qq;
	double anum,anum2;
	MyComplex temp, SLD, beta, rj;
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

	memset(pj, 0, (nlayers + 2) * sizeof(*pj));
	memset(SLDmatrix, 0, (nlayers + 2) * sizeof(*SLDmatrix));

	scale = coefP[1];
	bkg = fabs(coefP[4]);
	subrough = coefP[5];

	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 6;

	//fillout all the SLD's for all the layers
	for(ii = 1 ; ii < nlayers + 1 ; ii += 1){
		numtemp = 1.e-6 * ((100. - coefP[4*ii+4])/100.) * coefP[4*ii+3]+ (coefP[4*ii+4]*coefP[3]*1.e-6)/100.;		//sld of the layer
		*(SLDmatrix + ii) = 4 * PI * (numtemp  - (coefP[2] * 1e-6));
	}
	*(SLDmatrix) = 0;
	*(SLDmatrix + nlayers + 1) = 4 * PI * ((coefP[3] * 1e-6) - (coefP[2] * 1e-6));
	
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
		//set up an array for wavevectors
			try{
				SLDmatrixREP = new double [Vmullayers];
				pj_mul = new MyComplex [Vmullayers];
			} catch(...){
				err = NOMEM;
				goto done;
			}
			memset(pj_mul, 0, Vmullayers * sizeof(*pj_mul));
        
			for(ii=0; ii<Vmullayers;ii+=1){
				numtemp = (coefP[3]*1e-6*coefP[(4*ii)+offset+2]/100) +(1e-6 * ((100 - coefP[(4*ii)+offset+2])/100) * coefP[(4*ii)+offset+1]);		//sld of the layer
				*(SLDmatrixREP + ii) = 4 * PI * (numtemp  - (coefP[2] * 1e-6));
			}
	}
	
	for (j = 0; j < npoints ; j+=1) {
		//intialise the matrices
		memset(MRtotal, 0, sizeof(MRtotal));
		MRtotal[0][0].re = 1.;
		MRtotal[0][0].im = 0.;
		MRtotal[1][1].re = 1.;
		MRtotal[1][1].im = 0.;
		

		qq = xP[j] * xP[j] / 4;
		qq2.re = qq;
		qq2.im = 0;
		for(ii = 0; ii < nlayers+2 ; ii++)			//work out the wavevector in each of the layers
			pj[ii] = (*(SLDmatrix + ii) > qq) ? MyComplex(0, sqrt(fabs(qq - *(SLDmatrix + ii)))) : MyComplex(sqrt(qq - *(SLDmatrix + ii)), 0);
		
		//workout the wavevector in the toplayer of the multilayer, if it exists.
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
			memset(subtotal, 0, sizeof(subtotal));
			subtotal[0][0]=MyComplex(1,0);subtotal[1][1]=MyComplex(1,0);
			pj_mul[0] = (*(SLDmatrixREP)>qq) ? compsqrt(qq2-MyComplex(*SLDmatrixREP, 0)) : MyComplex(sqrt(qq-*SLDmatrixREP),0);
		}
		
		//now calculate reflectivities
		for(ii = 0 ; ii < nlayers+1 ; ii++){
			//work out the fresnel coefficients
			//this looks more complicated than it really is.
			//the reason it looks so convoluted is because if there is no complex part of the wavevector,
			//then it is faster to do the calc with real arithmetic then put it into a complex number.
			if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0 )
				rj = fres(pj[ii], pj_mul[0], coefP[offset+3]);
			else {
				if((pj[ii]).im == 0 && (pj[ii + 1]).im == 0){
					anum = (pj[ii]).re;
					anum2 = (pj[ii + 1]).re;
					rj.re = (ii == nlayers) ? 
					((anum - anum2) / (anum + anum2)) * exp(anum * anum2 * -2 * subrough * subrough)
					:
					((anum - anum2) / (anum + anum2)) * exp(anum * anum2 * -2 * coefP[4 * (ii + 1) + 5] * coefP[4 * (ii + 1) + 5]);
					rj.im = 0.;
				} else {
					rj = (ii == nlayers) ?
						((pj[ii] - pj[ii + 1])/(pj[ii] + pj[ii + 1])) * compexp(pj[ii] * pj[ii + 1] * -2 * subrough * subrough)
						:
						((pj[ii] - pj[ii + 1])/(pj[ii] + pj[ii + 1])) * compexp(pj[ii] * pj[ii + 1] * -2 * coefP[4 * (ii + 1) + 5] * coefP[4 * (ii + 1) + 5]);	
				};
			}

			//work out the beta for the (non-multi)layer
			temp.re = 0;
			temp.im = fabs(coefP[4 * ii + 2]);
			beta = (ii == 0)? oneC : compexp(pj[ii] * temp);

			//this is the characteristic matrix of a layer
			MI[0][0] = beta;
			MI[0][1] = rj * beta;
			MI[1][1] = oneC / beta;
			MI[1][0] = rj * MI[1][1];
			
			memcpy(temp2, MRtotal, sizeof(MRtotal));
						
			//multiply MR,MI to get the updated total matrix.			
			matmul(temp2, MI, MRtotal);

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

				memcpy(temp2, subtotal, sizeof(subtotal));

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

						memcpy(temp2, MRtotal, sizeof(MRtotal));
//						temp2[0][0] = MRtotal[0][0];
//						temp2[0][1] = MRtotal[0][1];
//						temp2[1][0] = MRtotal[1][0];
//						temp2[1][1] = MRtotal[1][1];

						matmul(temp2,MI,MRtotal);
					}
				} else {
					memcpy(temp2, MRtotal, sizeof(MRtotal));
//					temp2[0][0] = MRtotal[0][0];
//					temp2[0][1] = MRtotal[0][1];
//					temp2[1][0] = MRtotal[1][0];
//					temp2[1][1] = MRtotal[1][1];
					
					matmul(temp2,subtotal,MRtotal);
				};
			};
		};

		}
		
		den = compnorm(MRtotal[0][0]);
		num = compnorm(MRtotal[1][0]);
		answer=((num / den) * scale) + bkg;

		*yP++ = answer;
	}
	
done:
	if(pj != NULL)
		delete [] pj;
	if(pj_mul != NULL)
		delete[] pj_mul;
	if(SLDmatrix != NULL)
		delete[] SLDmatrix;
	if(SLDmatrixREP != NULL)
		delete[] SLDmatrixREP;

	return err;
}



 int
AbelesCalc_ImagAll(const double *coefP, double *yP, const double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep){
	int err = 0;
	int j;
	
	int ii=0,jj=0,kk=0;

	double scale,bkg,subrough;
	double num=0,den=0, answer=0;
	double anum, anum2;

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
		SLDmatrix = new MyComplex[nlayers + 2];
	} catch(...){
		err = NOMEM;
		goto done;
	}

	memset(pj, 0, (nlayers + 2) * sizeof(*pj));
	memset(SLDmatrix, 0, (nlayers + 2) * sizeof(*SLDmatrix));

	scale = coefP[1];
	bkg = coefP[6];
	subrough = coefP[7];
	sub= MyComplex(coefP[4]*1e-6, coefP[5]);
	super = MyComplex(coefP[2]*1e-6, coefP[3]);

	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 8;

	//fillout all the SLD's for all the layers
	for(ii=1; ii<nlayers+1;ii+=1)
		*(SLDmatrix + ii) = 4 * PI * (MyComplex(coefP[4 * ii + 5] * 1e-6, coefP[4 * ii + 6]) - super);

	*(SLDmatrix) = MyComplex(0,0);
	*(SLDmatrix + nlayers + 1) = 4 * PI * (sub - super);
	
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
	//set up an array for wavevectors
		try{
			SLDmatrixREP = new MyComplex[Vmullayers];
			pj_mul = new MyComplex[Vmullayers];
		} catch(...){
			err = NOMEM;
			goto done;
		}
		memset(pj_mul, 0, Vmullayers * sizeof(*pj_mul));
		memset(SLDmatrixREP, 0, Vmullayers * sizeof(*SLDmatrixREP));
		for(ii=0; ii<Vmullayers;ii+=1)
			*(SLDmatrixREP + ii) = 4 * PI * (MyComplex(coefP[(4 * ii) + offset + 1] * 1e-6, coefP[(4 * ii) + offset + 2])  - super);
	}


	for (j = 0; j < npoints; j++) {
		//intialise the matrices
		memset(MRtotal, 0, sizeof(MRtotal));
		MRtotal[0][0]=oneC;MRtotal[1][1]=oneC;

		qq2=MyComplex(xP[j]*xP[j]/4,0);

		for(ii=0; ii<nlayers+2 ; ii++){			//work out the wavevector in each of the layers
			pj[ii] = compsqrt(qq2-*(SLDmatrix+ii));
		}

		//workout the wavevector in the toplayer of the multilayer, if it exists.
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
			memset(subtotal, 0, sizeof(subtotal));
			subtotal[0][0]=MyComplex(1,0);subtotal[1][1]=MyComplex(1,0);
			pj_mul[0] = compsqrt(qq2-*SLDmatrixREP);
		}
		
		//now calculate reflectivities
		for(ii = 0 ; ii < nlayers+1 ; ii++){
			if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0 )
				rj = fres(pj[ii], pj_mul[0], coefP[offset+3]);
			else {
				if((pj[ii]).im == 0 && (pj[ii + 1]).im == 0){
					anum = (pj[ii]).re;
					anum2 = (pj[ii + 1]).re;
					rj.re = (ii == nlayers) ? 
					((anum - anum2) / (anum + anum2)) * exp(anum * anum2 * -2 * subrough * subrough)
					:
					((anum - anum2) / (anum + anum2)) * exp(anum * anum2 * -2 * coefP[4 * (ii + 1) + 7] * coefP[4 * (ii + 1) + 7]);
					rj.im = 0.;
				} else {
					rj = (ii == nlayers) ?
					((pj[ii] - pj[ii + 1])/(pj[ii] + pj[ii + 1])) * compexp(pj[ii] * pj[ii + 1] * -2 * subrough * subrough)
					:
					((pj[ii] - pj[ii + 1])/(pj[ii] + pj[ii + 1])) * compexp(pj[ii] * pj[ii + 1] * -2 * coefP[4 * (ii + 1) + 7] * coefP[4 * (ii + 1) + 7]);	
				};
			}
			

			//work out the beta for the (non-multi)layer
			beta = (ii==0)? oneC : compexp(pj[ii] * MyComplex(0,fabs(coefP[4*ii+4])));

			//this is the characteristic matrix of a layer
			MI[0][0]=beta;
			MI[0][1]=rj*beta;
			MI[1][1]=oneC/beta;
			MI[1][0]=rj*MI[1][1];

			memcpy(temp2, MRtotal, sizeof(MRtotal));

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
				((pj_mul[jj]-pj_mul[0])/(pj_mul[jj]+pj_mul[0]))* compexp((pj_mul[jj]*pj_mul[0])*-2*coefP[offset+3]*coefP[offset+3])
				:
				//otherwise it's the roughness of the layer below
				((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3]);
				
				
				//Beta's
				beta = compexp(MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]);

				MI[0][0]=beta;
				MI[0][1]=rj*beta;
				MI[1][1]=oneC/beta;
				MI[1][0]=rj*MI[1][1];

				memcpy(temp2, subtotal, sizeof(subtotal));

				matmul(temp2,MI,subtotal);
			};

			for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
				if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
					for(jj=0;jj<Vmullayers;jj++){
						beta = compexp((MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]));

						if(jj==Vmullayers-1){
							if(Vmulappend==nlayers){
								rj = ((pj_mul[Vmullayers-1]-pj[nlayers+1])/(pj_mul[Vmullayers-1]+pj[nlayers+1]))*compexp((pj_mul[Vmullayers-1]*pj[nlayers+1])*(-2*subrough*subrough));
							} else {
								rj = ((pj_mul[Vmullayers-1]-pj[Vmulappend+1])/(pj_mul[Vmullayers-1]+pj[Vmulappend+1]))* compexp((pj_mul[Vmullayers-1]*pj[Vmulappend+1])*(-2*coefP[4*(Vmulappend+1)+7]*coefP[4*(Vmulappend+1)+7]));
							};
						} else {
							rj = ((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3]);
						}
						
						MI[0][0]=beta;
						MI[0][1]=rj*beta;
						MI[1][1]=MyComplex(1,0)/MI[0][0];
						MI[1][0]=rj*MI[1][1];
						
						memcpy(temp2, MRtotal, sizeof(MRtotal));
//						temp2[0][0] = MRtotal[0][0];
//						temp2[0][1] = MRtotal[0][1];
//						temp2[1][0] = MRtotal[1][0];
//						temp2[1][1] = MRtotal[1][1];

						matmul(temp2,MI,MRtotal);
					}
				} else {
					memcpy(temp2, MRtotal, sizeof(MRtotal));
//					temp2[0][0] = MRtotal[0][0];
//					temp2[0][1] = MRtotal[0][1];
//					temp2[1][0] = MRtotal[1][0];
//					temp2[1][1] = MRtotal[1][1];
					
					matmul(temp2,subtotal,MRtotal);
				};
			};
		};

		}
		
		den= compnorm(MRtotal[0][0]);
		num=compnorm(MRtotal[1][0]);
		answer=(num/den);//(num*num)/(den*den);
		answer=(answer*scale)+fabs(bkg);

		*yP++ = answer;
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

	return err;
}

