/*
  FFTwrapper.c  -  A wrapper for Fast Fourier Transforms
  Author: Nasca O. Paul, Tg. Mures, Romania
  License: Public Domain
  
  This requires FFTW library (http://www.fftw.org) to compile and run.
*/

#include <math.h>
#include "FFTwrapper.h"

FFTwrapper::FFTwrapper(int fftsize_){
    fftsize=fftsize_;
    tmpfftdata1=new fftw_real[fftsize];
    tmpfftdata2=new fftw_real[fftsize];
#ifdef FFTW_VERSION_2
    planfftw=rfftw_create_plan(fftsize,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE|FFTW_IN_PLACE);
    planfftw_inv=rfftw_create_plan(fftsize,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE|FFTW_IN_PLACE);
#else
    planfftw=fftw_plan_r2r_1d(fftsize,tmpfftdata1,tmpfftdata1,FFTW_R2HC,FFTW_ESTIMATE);
    planfftw_inv=fftw_plan_r2r_1d(fftsize,tmpfftdata2,tmpfftdata2,FFTW_HC2R,FFTW_ESTIMATE);
#endif
};

FFTwrapper::~FFTwrapper(){
#ifdef FFTW_VERSION_2
    rfftw_destroy_plan(planfftw);
    rfftw_destroy_plan(planfftw_inv);
#else 
    fftw_destroy_plan(planfftw);
    fftw_destroy_plan(planfftw_inv);
#endif

    delete [] tmpfftdata1;
    delete [] tmpfftdata2;
};

/*
 * do the Fast Fourier Transform
 */
void FFTwrapper::smps2freqs(REALTYPE *smps,FFTFREQS freqs){
#ifdef FFTW_VERSION_2
    for (int i=0;i<fftsize;i++) tmpfftdata1[i]=smps[i];
    rfftw_one(planfftw,tmpfftdata1,tmpfftdata2);
    for (int i=0;i<fftsize/2;i++) {
	freqs.c[i]=tmpfftdata2[i];
	if (i!=0) freqs.s[i]=tmpfftdata2[fftsize-i];
    };
#else
    for (int i=0;i<fftsize;i++) tmpfftdata1[i]=smps[i];
    fftw_execute(planfftw);
    for (int i=0;i<fftsize/2;i++) {
	freqs.c[i]=tmpfftdata1[i];
	if (i!=0) freqs.s[i]=tmpfftdata1[fftsize-i];
    };
#endif
    tmpfftdata2[fftsize/2]=0.0;
};

/*
 * do the Inverse Fast Fourier Transform
 */
void FFTwrapper::freqs2smps(FFTFREQS freqs,REALTYPE *smps){
    tmpfftdata2[fftsize/2]=0.0;
#ifdef FFTW_VERSION_2
    for (int i=0;i<fftsize/2;i++) {
	tmpfftdata1[i]=freqs.c[i];
	if (i!=0) tmpfftdata1[fftsize-i]=freqs.s[i];
    };
    rfftw_one(planfftw_inv,tmpfftdata1,tmpfftdata2);
    for (int i=0;i<fftsize;i++) smps[i]=tmpfftdata2[i];
#else
    for (int i=0;i<fftsize/2;i++) {
	tmpfftdata2[i]=freqs.c[i];
	if (i!=0) tmpfftdata2[fftsize-i]=freqs.s[i];
    };
    fftw_execute(planfftw_inv);
    for (int i=0;i<fftsize;i++) smps[i]=tmpfftdata2[i];
#endif

};

void newFFTFREQS(FFTFREQS *f,int size){
    f->c=new REALTYPE[size];
    f->s=new REALTYPE[size];
    for (int i=0;i<size;i++){
	f->c[i]=0.0;f->s[i]=0.0;
    };
};
void deleteFFTFREQS(FFTFREQS *f){
    delete[] f->c;
    delete[] f->s;
    f->c=f->s=NULL;
};




