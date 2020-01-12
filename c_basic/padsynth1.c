/*
    Example implementation of the PADsynth basic algorithm
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
#include <assert.h>
#include <sys/time.h>

#include "../fftwrapper/FFTwrapper.h"
#define PI 3.14159265358979323846264338327950288


typedef struct {REALTYPE Re; REALTYPE Im;} complex;

/*
   fft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute fft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute fft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = -sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
void
fft2( complex *v, int n, complex *tmp )
{
	if (n > 1) {			/* otherwise, do nothing and return */
		int k, m;    complex z, w, *vo, *ve;
		ve = tmp; vo = tmp + n / 2;
		for (k = 0; k < n / 2; k++) {
			ve[k] = v[2 * k];
			vo[k] = v[2 * k + 1];
		}
		fft2( ve, n / 2, v );		/* FFT on even-indexed elements of v[] */
		fft2( vo, n / 2, v );		/* FFT on odd-indexed elements of v[] */
		for (m = 0; m < n / 2; m++) {
			w.Re = cos(2 * PI * m / (double)n);
			w.Im = -sin(2 * PI * m / (double)n);
			z.Re = w.Re * vo[m].Re - w.Im * vo[m].Im;	/* Re(w*vo[m]) */
			z.Im = w.Re * vo[m].Im + w.Im * vo[m].Re;	/* Im(w*vo[m]) */
			v[  m  ].Re = ve[m].Re + z.Re;
			v[  m  ].Im = ve[m].Im + z.Im;
			v[m + n / 2].Re = ve[m].Re - z.Re;
			v[m + n / 2].Im = ve[m].Im - z.Im;
		}
	}
	return;
}

/*
   ifft(v,N):
   [0] If N==1 then return.
   [1] For k = 0 to N/2-1, let ve[k] = v[2*k]
   [2] Compute ifft(ve, N/2);
   [3] For k = 0 to N/2-1, let vo[k] = v[2*k+1]
   [4] Compute ifft(vo, N/2);
   [5] For m = 0 to N/2-1, do [6] through [9]
   [6]   Let w.re = cos(2*PI*m/N)
   [7]   Let w.im = sin(2*PI*m/N)
   [8]   Let v[m] = ve[m] + w*vo[m]
   [9]   Let v[m+N/2] = ve[m] - w*vo[m]
 */
void
ifft2( complex *v, int n, complex *tmp )
{
	if (n > 1) {			/* otherwise, do nothing and return */
		int k, m;    complex z, w, *vo, *ve;
		ve = tmp; vo = tmp + n / 2;
		for (k = 0; k < n / 2; k++) {
			ve[k] = v[2 * k];
			vo[k] = v[2 * k + 1];
		}
		ifft2( ve, n / 2, v );		/* FFT on even-indexed elements of v[] */
		ifft2( vo, n / 2, v );		/* FFT on odd-indexed elements of v[] */
		for (m = 0; m < n / 2; m++) {
			w.Re = cos(2 * PI * m / (double)n);
			w.Im = sin(2 * PI * m / (double)n);
			z.Re = w.Re * vo[m].Re - w.Im * vo[m].Im;	/* Re(w*vo[m]) */
			z.Im = w.Re * vo[m].Im + w.Im * vo[m].Re;	/* Im(w*vo[m]) */
			v[  m  ].Re = ve[m].Re + z.Re;
			v[  m  ].Im = ve[m].Im + z.Im;
			v[m + n / 2].Re = ve[m].Re - z.Re;
			v[m + n / 2].Im = ve[m].Im - z.Im;
		}
	}
	return;
}




/* Random number generator */
REALTYPE RND() {
	return (rand() / (RAND_MAX + 1.0));
};

/* This is the profile of one harmonic
   In this case is a Gaussian distribution (e^(-x^2))
   The amplitude is divided by the bandwidth to ensure that the harmonic
   keeps the same amplitude regardless of the bandwidth */
REALTYPE profile(REALTYPE fi, REALTYPE bwi) {
	REALTYPE x = fi / bwi;
	x *= x;
	if (x > 14.71280603) return 0.0; //this avoids computing the e^(-x^2) where it's results are very close to zero
	return exp(-x) / bwi;
};

/*
    Inverse Fast Fourier Transform
    You may replace it with any IFFT routine
*/
void IFFT(int N, REALTYPE *freq_amp, REALTYPE *freq_phase, REALTYPE *smp) {
	FFTwrapper fft(N);
	FFTFREQS fftfreqs;
	newFFTFREQS(&fftfreqs, N / 2);

	for (int i = 0; i < N / 2; i++) {
		fftfreqs.c[i] = freq_amp[i] * cos(freq_phase[i]);
		fftfreqs.s[i] = freq_amp[i] * sin(freq_phase[i]);
	};
	fft.freqs2smps(fftfreqs, smp);
	deleteFFTFREQS(&fftfreqs);
};

/*
    Simple normalization function. It normalizes the sound to 1/sqrt(2)
*/
void normalize(int N, REALTYPE *smp) {
	int i;
	REALTYPE max = 0.0;
	for (i = 0; i < N; i++) if (fabs(smp[i]) > max) max = fabs(smp[i]);
	if (max < 1e-5) max = 1e-5;
	for (i = 0; i < N; i++) smp[i] /= max * 1.4142;
};


/*
    This is the implementation of PADsynth algorithm.
*/
void padsynth_basic_algorithm(
    int N, int samplerate, REALTYPE f, REALTYPE bw, int number_harmonics, REALTYPE *A, /*input data*/
    REALTYPE *smp /*output data*/)
{
	int i, nh;
	REALTYPE *freq_amp = new REALTYPE[N / 2];
	REALTYPE *freq_phase = new REALTYPE[N / 2];
	complex v[N / 2];
	complex smp2[N / 2];


	for (i = 0; i < N / 2; i++) freq_amp[i] = 0.0; //default, all the frequency amplitudes are zero

	for (nh = 1; nh < number_harmonics; nh++) { //for each harmonic
		REALTYPE bw_Hz;//bandwidth of the current harmonic measured in Hz
		REALTYPE bwi;
		REALTYPE fi;
		bw_Hz = (pow(2.0, bw / 1200.0) - 1.0) * f * nh;

		bwi = bw_Hz / (2.0 * samplerate);
		fi = f * nh / samplerate;
		for (i = 0; i < N / 2; i++) {
			REALTYPE hprofile;
			hprofile = profile((i / (REALTYPE)N) - fi, bwi);
			freq_amp[i] += hprofile * A[nh];
			v[i].Im = freq_amp[i];
		};
	};

	//Add random phases
	srandom(time(0));
	for (i = 0; i < N / 2; i++) {
		freq_phase[i] = RND() * 2.0 * PI;
		v[i].Re = freq_phase[i];
		smp2[i].Re = 0;
		smp2[i].Im = 0;
	};

	IFFT(N, freq_amp, freq_phase, smp);
	// ifft2(v,N/2,smp2);
	// for (i=0;i<N/2;i++) {
	// 	// printf("%2.3f ",smp[i]);
	// 	// printf("%2.3f\n",sqrtf(pow(v[i].Re,2)+pow(v[i].Im,2)));
	// 	smp[i] = sqrtf(pow(v[i].Re,2)+pow(v[i].Im,2));
	// }
	normalize(N, smp);

	delete [] freq_amp;
	delete [] freq_phase;
};

void sleep_ms(int milliseconds) // cross-platform sleep function
{
#ifdef WIN32
	Sleep(milliseconds);
#elif _POSIX_C_SOURCE >= 199309L
	struct timespec ts;
	ts.tv_sec = milliseconds / 1000;
	ts.tv_nsec = (milliseconds % 1000) * 1000000;
	nanosleep(&ts, NULL);
#else
	usleep(milliseconds * 1000);
#endif
}

float timedifference_msec(struct timeval t0, struct timeval t1)
{
	return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}


#define N (800000)

int main() {
	// TODO try to generate only the number of samples to keep it in phase

	srandom(time(0));
	REALTYPE A[128]; A[0] = 0.0; //A[0] is not used
	REALTYPE A2[128]; A2[0] = 0.0; //A[0] is not used
	REALTYPE num;
	num = 0;
	REALTYPE Nstop;
	REALTYPE Nstopbreak;
	short int value;
	REALTYPE sample[N];
	short int isample[N];
	struct timeval t0;
	struct timeval t1;
	float elapsed;
	FILE *myfile;
	double myvariable;

	// initialize


	// 100 millisecond intervals
	REALTYPE freqs[] = {130.81,  164.81, 220, 246.94};
	int number_harmonics = 128;
	int waviness = 10;
	myfile = fopen("values.txt", "r");
	for (int i = 0; i < 4; i++) {
		fscanf(myfile, "%f", &freqs[i]);
	}
	fscanf(myfile, "%d", &number_harmonics);
	fscanf(myfile, "%d", &waviness);
	fclose(myfile);


	REALTYPE lowfreq = 10000;
	REALTYPE highfreq = 0;
	REALTYPE numfreq = 0;
	for (int j1 = 0; j1 < 1000; j1++) {
		gettimeofday(&t0, 0);

		for (int i = 1; i < 128; i++) {
			A[i] = 1.0 / i;
			if ((i % 2) == 0) A[i] *= 2.0;
		};
		numfreq = 0;
		highfreq = 0;
		lowfreq = 10000;
		for (int i = 0; i < 4; i++) {
			if (freqs[i] > 0) {
				numfreq++;
			}
			if (freqs[i] > 0 && freqs[i] < lowfreq) {
				lowfreq = freqs[i];
			}
			if (freqs[i] > 0 && freqs[i] > highfreq) {
				highfreq = freqs[i];
			}
		}
		fprintf(stderr, "numfreq = %f\n", numfreq);
		Nstop = 44100.0 / highfreq * 500;
		while (Nstop > N / 2) {
			Nstop = Nstop - 44100.0 / highfreq ;
		}

		for (int i = 0; i < N; i++) {
			sample[i] = 0;
		}

		fprintf(stderr, "Nstop = %f\n", Nstop);
		padsynth_basic_algorithm(Nstop, 44100, freqs[0], waviness, number_harmonics, A, sample);
		for (int i = 0; i < Nstop; i++) {
			// isample[i] = (int)(sample[i] * 32768.0/2 + sample2[i]* 32768.0/2);
			isample[i] = (int)(sample[i] * 32768.0 / numfreq);
		}
		for (int k = 1; k < 4; k++) {
			if (freqs[k] > 0) {
				padsynth_basic_algorithm(Nstop, 44100, freqs[k], waviness, number_harmonics, A, sample);
				for (int i = 0; i < Nstop; i++) {
					// isample[i] = (int)(sample[i] * 32768.0/2 + sample2[i]* 32768.0/2);
					isample[i] = isample[i] + (int)(sample[i] * 32768.0 / numfreq);
				}
			}
		}

		gettimeofday(&t1, 0);
		elapsed = timedifference_msec(t0, t1);
		fprintf(stderr, "Code executed in %f milliseconds.\n", elapsed);

		// fwrite(isample,Nstop,2,stdout);
		Nstopbreak = 44100.0 / highfreq;
		for (int i = 0; i < Nstop; i++) {
			fwrite(&isample[i], 1, 2, stdout);
			// // check at each Nstop to see if we have new info
			if (i > Nstopbreak) {
				REALTYPE currentvalue = freqs[0] + freqs[1] + freqs[2] + freqs[3] + number_harmonics + waviness;
				myfile = fopen("values.txt", "r");
				for (int i = 0; i < 4; i++) {
					fscanf(myfile, "%f", &freqs[i]);
				}
				fscanf(myfile, "%d", &number_harmonics);
				fscanf(myfile, "%d", &waviness);
				fclose(myfile);
				REALTYPE nextvalue = freqs[0] + freqs[1] + freqs[2] + freqs[3] + number_harmonics + waviness;
				Nstopbreak += 44100.0 / highfreq;
				if (nextvalue != currentvalue) {
					fprintf(stderr, "breaking\n");
					break;
				}
			}
		}
		// sleep_ms( Nstop/44100 * 1000-250);

		fprintf(stderr, "j: %d\n", j1);
		int sleeptime =  Nstop / 44100 * 1000 - j1 - 200;
		fprintf(stderr, "Sleeping %d ms\n", sleeptime);
		// sleep_ms(sleeptime);
	}



//    for (int i=0;i<N;i++) {
	// 	fwrite(&isample[N-i-1],1,2,stdout);
	// }
//    for (int i=0;i<N-2;i+=2) {
//    	for (int j=0;j<2;j++) {
//    		subset[j] = isample[(N-i+j)-1];
//    	}
	// 	fwrite(subset,2,2,stdout);
	// }
	// FILE *f=fopen("sample.raw","w");fwrite(isample,N,2,f);fclose(f);
};


