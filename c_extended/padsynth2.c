/*
    Example implementation of the PADsynth extended algorithm
    By: Nasca O. Paul, Tg. Mures, Romania
    This implementation and the algorithm are released under Public Domain
    Feel free to use it into your projects or your products ;-)

    This implementation is tested under GCC/Linux, but it's 
    very easy to port to other compiler/OS.
    
    P.S.
    Please note, that IFFT function depends on the FFTW library, so if you want 
    to use into commercial products, you must replace it with your IFFT routine
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../fftwrapper/FFTwrapper.h"
#define PI 3.14159265358979;

/* Random number generator */
REALTYPE RND(){
    return (rand()/(RAND_MAX+1.0));
};

/* This is the profile of one harmonic
   In this case is a Gaussian distribution (e^(-x^2))
   The amplitude is divided by the bandwidth to ensure that the harmonic
   keeps the same amplitude regardless of the bandwidth */
REALTYPE profile(REALTYPE fi,REALTYPE bwi){
    REALTYPE x=fi/bwi;
    x*=x;
    if (x>14.71280603) return 0.0;//this avoids computing the e^(-x^2) where it's results are very close to zero
    return exp(-x)/bwi;
};

/*
    Inverse Fast Fourier Transform
    You may replace it with any IFFT routine
*/
void IFFT(int N,REALTYPE *freq_amp,REALTYPE *freq_phase,REALTYPE *smp){
    FFTwrapper fft(N);
    FFTFREQS fftfreqs;
    newFFTFREQS(&fftfreqs,N/2);

    for (int i=0;i<N/2;i++){
	fftfreqs.c[i]=freq_amp[i]*cos(freq_phase[i]);
	fftfreqs.s[i]=freq_amp[i]*sin(freq_phase[i]);
    };
    fft.freqs2smps(fftfreqs,smp);
    deleteFFTFREQS(&fftfreqs);
};

/*
    Simple normalization function. It normalizes the sound to 1/sqrt(2)
*/
void normalize(int N,REALTYPE *smp){
    int i;
    REALTYPE max=0.0;
    for (i=0;i<N;i++) if (fabs(smp[i])>max) max=fabs(smp[i]);
    if (max<1e-5) max=1e-5;
    for (i=0;i<N;i++) smp[i]/=max*1.4142;
};

/*
    The relF function returns the relative frequency of the N'th harmonic
    to the fundamental frequency.
*/
REALTYPE relF(int N){
    return (N*(1.0+N*0.1));
};

/*
    This is the implementation of PADsynth algorithm.
*/
void padsynth_extended_algorithm(
	/*input data*/
	int N,
	int samplerate,
	REALTYPE f,
	REALTYPE bw,
	REALTYPE bwscale,
	int number_harmonics,
	REALTYPE *A, 
	/*output data*/
        REALTYPE *smp )
{
    int i,nh;
    REALTYPE *freq_amp=new REALTYPE[N/2];
    REALTYPE *freq_phase=new REALTYPE[N/2];
    
    for (i=0;i<N/2;i++) freq_amp[i]=0.0;//default, all the frequency amplitudes are zero

    for (nh=1;nh<number_harmonics;nh++){//for each harmonic
	REALTYPE bw_Hz;//bandwidth of the current harmonic measured in Hz
        REALTYPE bwi;
	REALTYPE fi;

        bw_Hz=(pow(2.0,bw/1200.0)-1.0)*f*pow(relF(nh),bwscale);
	
	bwi=bw_Hz/(2.0*samplerate);
	fi=f*relF(nh)/samplerate;
	for (i=0;i<N/2;i++){
	    REALTYPE hprofile;
	    hprofile=profile((i/(REALTYPE)N)-fi,bwi);
	    freq_amp[i]+=hprofile*A[nh];
	};
    };
    
    //Add random phases
    for (i=0;i<N/2;i++){
	freq_phase[i]=RND()*2.0*PI;
    };
    
    IFFT(N,freq_amp,freq_phase,smp);
    normalize(N,smp);
    
    delete [] freq_amp;
    delete [] freq_phase;
};

#define N (262144)
REALTYPE sample[N];

#define number_harmonics 128
int main(){
    srandom(time(0));
    REALTYPE A[number_harmonics];A[0]=0.0;//A[0] is not used
    for (int i=1;i<number_harmonics;i++) A[i]=1.0/i;
    padsynth_extended_algorithm(N,44100,101.0,10.0,10.9,number_harmonics,A,sample);

    /* Output the data to the 16 bit, mono raw file */
    short int isample[N];
    for (int i=0;i<N;i++) isample[i]=(int)(sample[i]*32768.0);
    FILE *f=fopen("sample.raw","w");fwrite(isample,N,2,f);fclose(f);
};


