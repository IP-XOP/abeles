/*
 *  Threadcalculator.cpp
 *  Abeles
 *
 Copyright Andrew Nelson and ANSTO 2009
 *
 */
#include "XOPStandardHeaders.h"
#include "RefCalculator.h"
#include <complex>
#include "Abeles.h"

#ifdef MACIGOR
	#include <pthread.h>
#endif

#ifdef WINIGOR
#include "pthread.h"
#include "sched.h"
#include "semaphore.h"
#endif

using namespace std;


void *AbelesThreadWorker(void *arg){
	int err = NULL;
	refCalcParm *p = (refCalcParm *) arg;
	err = AbelesCalcAll(p->coefP, p->yP, p->xP, p->npoints, p->Vmullayers, p->Vappendlayer, p->Vmulrep);
	
	pthread_exit((void*) &err);
	return NULL;
}


void *AbelesImagThreadWorker(void *arg){
	int err = NULL;
	refCalcParm *p = (refCalcParm *) arg;
	err = AbelesCalc_ImagAll(p->coefP, p->yP, p->xP, p->npoints, p->Vmullayers, p->Vappendlayer, p->Vmulrep);
	
	pthread_exit((void*) &err);
	return NULL;
}


int 
AbelesCalcAll(const double *coefP, double *yP, const double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep){
	int err = 0;
	int j;
	
	int ii=0,jj=0,kk=0;

	double scale,bkg,subrough;
	double num=0, den=0, answer=0;
	complex<double> temp, SLD, beta, rj;
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
	bkg = coefP[4];
	subrough = coefP[5];

	// offset tells us where the multilayers start.
	offset = 4 * nlayers + 6;

	// fillout all the SLD's for all the layers
	for(ii = 1 ; ii < nlayers + 1 ; ii += 1){
		numtemp = 1.e-6 * ((100. - coefP[4*ii+4])/100.) * coefP[4*ii+3]+ (coefP[4*ii+4]*coefP[3]*1.e-6)/100.;		//sld of the layer
		*(SLDmatrix + ii) = 4 * PI * (numtemp  - (coefP[2] * 1e-6));
	}
	*(SLDmatrix) = 0;
	*(SLDmatrix + nlayers + 1) = 4 * PI * ((coefP[3] * 1e-6) - (coefP[2] * 1e-6));
	
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
		// set up an array for wavevectors
			try{
				SLDmatrixREP = new double [Vmullayers];
				pj_mul = new complex<double> [Vmullayers];
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
		// intialise the matrices
		memset(MRtotal, 0, sizeof(MRtotal));
        MRtotal[0][0] = oneC;
        MRtotal[1][1] = oneC;

		qq2 = complex<double>(xP[j] * xP[j] / 4, 0);
		for(ii = 0; ii < nlayers + 2 ; ii++)
            // work out the wavevector in each of the layers
            pj[ii] = std::sqrt(qq2 - *(SLDmatrix + ii));
		
		// workout the wavevector in the toplayer of the multilayer, if it exists.
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >=0){
			memset(subtotal, 0, sizeof(subtotal));
			subtotal[0][0] = complex<double>(1, 0);
            subtotal[1][1] = complex<double>(1, 0);
            pj_mul[0] = std::sqrt(qq2 - *SLDmatrixREP);
		}
		
		// now calculate reflectivities
		for(ii = 0 ; ii < nlayers + 1 ; ii++){
			// work out the fresnel coefficients
			// this looks more complicated than it really is.
			// the reason it looks so convoluted is because if there is no complex part of the wavevector,
			// then it is faster to do the calc with real arithmetic then put it into a complex number.
			if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0 )
				rj = fres(pj[ii], pj_mul[0], coefP[offset+3]);
			else {
                rj = (ii == nlayers) ?
                fres(pj[ii], pj[ii + 1], subrough)
                    :
                fres(pj[ii], pj[ii + 1], coefP[4 * (ii + 1) + 5]);
			}

			//work out the beta for the (non-multi)layer
            temp = complex<double>(0, fabs(coefP[4 * ii + 2]));
            beta = (ii == 0)? oneC : std::exp(pj[ii] * temp);

			//this is the characteristic matrix of a layer
			MI[0][0] = beta;
			MI[0][1] = rj * beta;
			MI[1][1] = oneC / beta;
			MI[1][0] = rj * MI[1][1];
			
			memcpy(temp2, MRtotal, sizeof(MRtotal));
						
			//multiply MR,MI to get the updated total matrix.			
			matmul(temp2, MI, MRtotal);

		if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0){
		    // workout the wavevectors in each of the layers
			for(jj=1 ; jj < Vmullayers; jj++){
                pj_mul[jj] = std::sqrt(qq2 - *(SLDmatrixREP+jj));
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
				MI[1][1] = oneC/beta;
				MI[1][0] = rj * MI[1][1];

				memcpy(temp2, subtotal, sizeof(subtotal));

				matmul(temp2,MI,subtotal);
			};

			for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
				if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
					for(jj=0;jj<Vmullayers;jj++){
                        beta = std::exp((complex<double>(0, fabs(coefP[4*jj+offset])) * pj_mul[jj]));

						if(jj == Vmullayers-1){
							if(Vmulappend == nlayers){
                                rj = fres(pj_mul[Vmullayers-1], pj[nlayers+1], subrough);
							} else {
                                rj = fres(pj_mul[Vmullayers-1], pj[Vmulappend+1], coefP[4*(Vmulappend+1)+5]);
							};
						} else {
                            rj = fres(pj_mul[jj], pj_mul[jj+1], coefP[4*(jj+1)+offset+3]);
						}

						MI[0][0] = beta;
						MI[0][1] = (rj * beta);
						MI[1][1] = oneC / beta;
						MI[1][0] = (rj * MI[1][1]);

						memcpy(temp2, MRtotal, sizeof(MRtotal));
                        matmul(temp2,MI,MRtotal);
					}
				} else {
					memcpy(temp2, MRtotal, sizeof(MRtotal));
					matmul(temp2,subtotal,MRtotal);
				};
			};
		};

		}
		
        num = std::norm(MRtotal[1][0]);
        den = std::norm(MRtotal[0][0]);
        answer = (num / den);
        answer = (answer * scale) + bkg;
        
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
AbelesCalc_ImagAll(const double *coefP, double *yP, const double *xP, long npoints, int Vmullayers, int Vmulappend, int Vmulrep){
	int err = 0;
	int j;
	
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
		SLDmatrix = new complex<double>[nlayers + 2];
	} catch(...){
		err = NOMEM;
		goto done;
	}

	memset(pj, 0, (nlayers + 2) * sizeof(*pj));
	memset(SLDmatrix, 0, (nlayers + 2) * sizeof(*SLDmatrix));

	scale = coefP[1];
	bkg = coefP[6];
	subrough = coefP[7];
	sub = complex<double>(coefP[4] * 1e-6, coefP[5] + 1.e-100);
	super = complex<double>(coefP[2] * 1e-6, 0);

	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 8;

	//fillout all the SLD's for all the layers
	for(ii = 1; ii < nlayers + 1; ii += 1)
		SLDmatrix[ii] = 4 * PI * (complex<double>(coefP[4 * ii + 5] * 1.e-6, coefP[4 * ii + 6] + 1.e-100) - super);

	*(SLDmatrix) = complex<double>(0, 0);
	*(SLDmatrix + nlayers + 1) = 4 * PI * (sub - super);
	
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
		for(ii = 0; ii < Vmullayers; ii += 1)
			*(SLDmatrixREP + ii) = 4 * PI * (complex<double>(coefP[(4 * ii) + offset + 1] * 1e-6, coefP[(4 * ii) + offset + 2] + 1.e-100)  - super);
	}

	for (j = 0; j < npoints; j++) {
		//intialise the matrices
		memset(MRtotal, 0, sizeof(MRtotal));
		MRtotal[0][0] = oneC;
        MRtotal[1][1] = oneC;

		qq2 = complex<double>(xP[j]*xP[j]/4, 0);

		for(ii = 0; ii < nlayers + 2 ; ii++){			//work out the wavevector in each of the layers
            pj[ii] = std::sqrt(qq2 - *(SLDmatrix + ii));
		}

		//workout the wavevector in the toplayer of the multilayer, if it exists.
		if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
			memset(subtotal, 0, sizeof(subtotal));
			subtotal[0][0] = complex<double>(1, 0);
            subtotal[1][1] = complex<double>(1, 0);
            pj_mul[0] = std::sqrt(qq2 - *SLDmatrixREP);
		}
		
		//now calculate reflectivities
		for(ii = 0 ; ii < nlayers+1 ; ii++){
            if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0 ){
				rj = fres(pj[ii], pj_mul[0], coefP[offset+3]);
            } else {
                rj = (ii == nlayers) ?
                fres(pj[ii], pj[ii + 1], subrough)
                :
                fres(pj[ii], pj[ii + 1], coefP[4 * (ii + 1) + 7]);
            };

            //work out the beta for the (non-multi)layer
            beta = (ii == 0)? oneC : std::exp(pj[ii] * complex<double>(0, fabs(coefP[4 * ii + 4])));

            //this is the characteristic matrix of a layer
            MI[0][0] = beta;
            MI[0][1] = rj * beta;
            MI[1][1] = oneC / beta;
            MI[1][0] = rj * MI[1][1];

            memcpy(temp2, MRtotal, sizeof(MRtotal));

            // multiply MR, MI to get the updated total matrix.
            matmul(temp2, MI, MRtotal);

            if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0){
                // workout the wavevectors in each of the layers
                for(jj=1 ; jj < Vmullayers; jj++)
                    pj_mul[jj] = std::sqrt(qq2 - *(SLDmatrixREP + jj));

                //work out the fresnel coefficients
                for(jj = 0 ; jj < Vmullayers; jj++){
                    rj = (jj == Vmullayers - 1) ?
                    //if you're in the last layer then the roughness is the roughness of the top
                    fres(pj_mul[jj], pj_mul[0], coefP[offset+3])
                    :
                    //otherwise it's the roughness of the layer below
                    fres(pj_mul[jj],  pj_mul[jj+1], coefP[4*(jj+1)+offset+3]);
                    
                    // Beta's
                    beta = std::exp(complex<double>(0, fabs(coefP[4 * jj + offset])) * pj_mul[jj]);

                    MI[0][0] = beta;
                    MI[0][1] = rj * beta;
                    MI[1][1] = oneC / beta;
                    MI[1][0] = rj * MI[1][1];

                    memcpy(temp2, subtotal, sizeof(subtotal));

                    matmul(temp2,MI,subtotal);
                };

                for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
                    if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
                        for(jj = 0; jj < Vmullayers; jj++){
                            beta = std::exp((complex<double>(0, fabs(coefP[4*jj+offset])) * pj_mul[jj]));

                            if(jj == Vmullayers - 1){
                                if(Vmulappend == nlayers){
                                    rj = fres(pj_mul[Vmullayers-1], pj[nlayers+1], subrough);
                                } else {
                                    rj = fres(pj_mul[Vmullayers-1], pj[Vmulappend+1], coefP[4*(Vmulappend+1) + 7]);
                                };
                            } else {
                                rj = fres(pj_mul[jj], pj_mul[jj+1], coefP[4*(jj+1)+offset+3]);
                            }
                            
                            MI[0][0] = beta;
                            MI[0][1] = rj * beta;
                            MI[1][1] = complex<double>(1, 0) / MI[0][0];
                            MI[1][0] = rj * MI[1][1];
                            
                            memcpy(temp2, MRtotal, sizeof(MRtotal));
                            matmul(temp2,MI,MRtotal);
                        }
                    } else {
                        memcpy(temp2, MRtotal, sizeof(MRtotal));
                        matmul(temp2,subtotal,MRtotal);
                    };
                };
            };
        };
        num = std::norm(MRtotal[1][0]);
        den = std::norm(MRtotal[0][0]);
        answer = (num / den);
        answer = (answer * scale) + bkg;

        *yP++ = answer;
    };
	
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
