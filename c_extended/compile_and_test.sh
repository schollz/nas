rm -f a.out ; g++ *.c ../fftwrapper/*.cpp -lm -lfftw3 ; ./a.out ; 
sox -r 44100 -sw -c 1 sample.raw sample.wav 
#play sample.wav
rm -f sample.raw 

