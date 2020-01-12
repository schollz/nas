build:
	g++ *.c ../fftwrapper/*.cpp -lm -lfftw3 -lpthread -o pad1

play:
	./pad1 | ../alsa/alsastream 44100 1 1