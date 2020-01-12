rm -f a.out ; g++ *.c ../fftwrapper/*.cpp -lm -lfftw3 ; ./a.out ; 
for file in note??
do
    echo Converting $file 
    sox -r 44100 -sw -c 1 -t raw $file $file.wav 
done
rm -f note??

