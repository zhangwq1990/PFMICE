mpif90 -r8 -fpconstant -O3 -132 -vec-report0 -c fft.f
ar qc libfft.a fft.o
cp libfft.a lib
