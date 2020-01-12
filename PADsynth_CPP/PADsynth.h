/*
    PADsynth implementation as ready-to-use C++ class.
    By: Nasca O. Paul, Tg. Mures, Romania
    This implementation and the algorithm are released under Public Domain
    Feel free to use it into your projects or your products ;-)

    This implementation is tested under GCC/Linux, but it's 
    very easy to port to other compiler/OS.

    Please notice that you have to implement the virtual method IFFT() with your IFFT routine.
    
*/

#ifndef PADSYNTH_H
#define PADSYNTH_H

#ifndef REALTYPE
#define REALTYPE float
#endif

class PADsynth{
    public:
	/*  PADsynth:
		N                - is the samplesize (eg: 262144)
		samplerate 	 - samplerate (eg. 44100)
		number_harmonics - the number of harmonics that are computed */
	PADsynth(int N_,int samplerate_,int number_harmonics_);

	~PADsynth();

	/* set the amplitude of the n'th harmonic */
	void setharmonic(int n,REALTYPE value);

	/* get the amplitude of the n'th harmonic */
	REALTYPE getharmonic(int n);

	/*  synth() generates the wavetable
	    f		- the fundamental frequency (eg. 440 Hz)
	    bw		- bandwidth in cents of the fundamental frequency (eg. 25 cents)
	    bwscale	- how the bandwidth increase on the higher harmonics (recomanded value: 1.0)
	    *smp	- a pointer to allocated memory that can hold N samples */
	void synth(REALTYPE f,REALTYPE bw,
	    REALTYPE bwscale,REALTYPE *smp);
    protected:
	int N;			//Size of the sample

	/* IFFT() - inverse fast fourier transform
	   YOU MUST IMPLEMENT THIS METHOD!
	   *freq_real and *freq_imaginary represents the real and the imaginary part of the spectrum, 
	   The result should be in *smp array.
	   The size of the *smp array is N and the size of the freq_real and freq_imaginary is N/2 */
	virtual void IFFT(REALTYPE *freq_real,REALTYPE *freq_imaginary,REALTYPE *smp)=0;


	/* relF():
	    This method returns the N'th overtone's position relative 
	    to the fundamental frequency.
	    By default it returns N.
	    You may override it to make metallic sounds or other 
	    instruments where the overtones are not harmonic.  */
	virtual REALTYPE relF(int N);

	/* profile():
	    This is the profile of one harmonic
	    In this case is a Gaussian distribution (e^(-x^2))
            The amplitude is divided by the bandwidth to ensure that the harmonic
	    keeps the same amplitude regardless of the bandwidth */
	virtual REALTYPE profile(REALTYPE fi, REALTYPE bwi);

	/* RND() - a random number generator that 
	    returns values between 0 and 1
	*/
	virtual REALTYPE RND();

    private:
	REALTYPE *A;		//Amplitude of the harmonics
	REALTYPE *freq_amp;	//Amplitude spectrum
	int samplerate;
	int number_harmonics;
};

#endif

