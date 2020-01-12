#ifndef FFT_WRAPPER_H
#define FFT_WRAPPER_H

#define REALTYPE float

struct FFTFREQS{
    REALTYPE *s,*c;//sine and cosine components
};

#ifdef FFTW_VERSION_2

#include <fftw.h>

/* If you got error messages about rfftw.h, replace the next include  line with "#include <srfftw.h>"
or with "#include <drfftw.h> (if one doesn't work try the other). It may be necessary to replace
the <fftw.h> with <dfftw.h> or <sfftw.h>. If the neither one doesn't work, 
please install latest version of fftw(recomanded from the sources) from www.fftw.org.
If you'll install fftw3 you need to change the Makefile.inc
Hope all goes right." */
#include <rfftw.h>

#else

#include <fftw3.h>
#define fftw_real double
#define rfftw_plan fftw_plan
#endif

void newFFTFREQS(FFTFREQS *f,int size);
void deleteFFTFREQS(FFTFREQS *f);
class FFTwrapper{
    public:
	FFTwrapper(int fftsize_);
	~FFTwrapper();
	void smps2freqs(REALTYPE *smps,FFTFREQS freqs);
	void freqs2smps(FFTFREQS freqs,REALTYPE *smps);
    private:
	int fftsize;
	fftw_real *tmpfftdata1,*tmpfftdata2;
	rfftw_plan planfftw,planfftw_inv;
};
#endif

