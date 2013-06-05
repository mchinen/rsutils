/****************************************************************
 
    rs_fft.h

    Complex and Real FFT/IFFT

    Michael Chinen
    http://roughsoft.com

    Copyright 2012 Roughsoft LLC
    License: MIT

****************************************************************/

#ifndef __RS_FFT_H__
#define __RS_FFT_H__


#ifdef __cplusplus
extern "C" {
#endif


/// This function is used to do transformations to and from purely
/// real values (e.g. audio data).
/// When used as an FFT, this function takes N * 2 real values 
/// and returns the first half plus the nyquist bin (N + 1) of its
/// complex spectrum in-place.  x[1] contains the N+1th real value.
/// When used as an IFFT, takes N complex values that represent
/// the first half of the spectrum.  x[1] should contain the N+1th bin
/// @param N should be set so that there are 2 * N real values in x
/// @param forward if 1, does an fft.  If 0, does an ifft.
void rs_rfft(double* x, int N, int forward);


/// A radix-2 decimation in time FFT that takes a complex input
/// @param x an array with N interleaving real and imaginary values
/// @param N the number of complex values in x
/// @param forward if 1 does an fft.  If 0, does an ifft.
void rs_cfft(double* x, int N, int forward);


#ifdef __cplusplus
};
#endif

#endif
