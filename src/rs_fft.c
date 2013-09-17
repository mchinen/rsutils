/****************************************************************
 
    rs_fft.h

    Complex and Real FFT/IFFT

    Michael Chinen
    http://roughsoft.com

    Copyright 2012 Roughsoft LLC
    License: MIT

****************************************************************/


//get M_PI
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

// x has N * 2 values (moved in pairs.
static void rs_reverse(double* x, int N)
{
   double tempi, tempr;
   unsigned int i;   
   static int lastN = -1;
   static unsigned int* mirror;

   // cache a mirror index since doing this on the fly will
   // yeild O(N logN) just for the the reverse.
   if (lastN != N) {
      unsigned int j;
      int log2N = 0;
      unsigned int temp = N;
      if (mirror)
         free(mirror);

      // malloc and memset to zero
      mirror = (unsigned int*) calloc(sizeof(unsigned int) * N, 1);

      while (temp >>= 1) {
         log2N++;
      }
      for (i = 0; i < N; i++) {
         temp = i;
         for (j = 0; j < log2N; j++) {
            // we can probably make this more efficient if needed, but since this is cached..
            mirror[i] = (mirror[i] << 1) | (temp & 1);
            temp = temp >> 1;
         }
      }

      lastN = N;
   }

   for (i = 0; i < N; i++) {
      if (i < mirror[i]) {
         tempr = x[i * 2];
         tempi = x[i * 2 + 1];
         x[i * 2]     = x[mirror[i] * 2];
         x[i * 2 + 1] = x[mirror[i] * 2 + 1];
         x[mirror[i] * 2]     = tempr;
         x[mirror[i] * 2 + 1] = tempi;
      }
   }
}

static double* fetch_twiddle_array(int N, int forward)
{
   int i;
   double** twiddleHandle;
   int* lastNPointer;

   // inverse and forward use seperate twiddle values, so we use handles to manage this
   static double* twiddleForward;
   static double* twiddleInverse;
   static int lastNForward = -1;
   static int lastNInverse = -1;

   twiddleHandle = forward ? &twiddleForward : &twiddleInverse;
   lastNPointer  = forward ? &lastNForward   : &lastNInverse;

   if (*lastNPointer != N) {
      if (*twiddleHandle)
         free(*twiddleHandle);
     
      *twiddleHandle = malloc(sizeof(double) * N);
      for (i = 0; i < N  / 2; i++) {
         // we could use a recursive solution that multiplies against e^(pi/N)
         // but I suspect that is lossy.
         // since we are caching it, do it at high accuracy once.
         // note that we are computing a phi of -j2pik/N.
         // we only do half because we can use symmetry to calculate the rest.
         // also note the negative sign in front of j means we walk backwards,
         // but the ifft walks forward, which is why we have two arrays.
         // symmetry could be used there too, but it's not as simple for the client.
         // TODO: see if the -2 * sin^2 phi substitute for cos on grounds
         // of lossiness has grounds here.  I suspect it doesn't for our N.
         (*twiddleHandle)[i * 2]     = cos(M_PI * (forward ? -i : i) / (N / 2));
         (*twiddleHandle)[i * 2 + 1] = sin(M_PI * (forward ? -i : i) / (N / 2));
      }
      *lastNPointer = N;
   }
   return *twiddleHandle;
}

// this is kind of silly but the rfft uses a larger N than the cfft.
// so we use a duplicate function so the two static arrays can coexist.
// the only difference functionally is that the divisor for theta in sin/cos
// is N instead of N / 2.
// TODO: take advantange of the fact that cfft's twiddles are a subset of rfft's
// if the spacing is 2
static double* fetch_twiddle_array_rfft(int N, int forward)
{
   int i;
   double** twiddleHandle;
   int* lastNPointer;

   // inverse and forward use seperate twiddle values, so we use handles to manage this
   static double* twiddleForward;
   static double* twiddleInverse;
   static int lastNForward = -1;
   static int lastNInverse = -1;

   twiddleHandle = forward ? &twiddleForward : &twiddleInverse;
   lastNPointer  = forward ? &lastNForward   : &lastNInverse;

   if (*lastNPointer != N) {
      if (*twiddleHandle)
         free(*twiddleHandle);
     
      *twiddleHandle = malloc(sizeof(double) * N);
      for (i = 0; i < N  / 2; i++) {
         // we could use a recursive solution that multiplies against e^(pi/N)
         // but I suspect that is lossy.
         // since we are caching it, do it at high accuracy once.
         // TODO: see if the -2 * sin^2 phi substitute for cos on grounds
         // of lossiness has grounds here.  I suspect it doesn't for our N.
         (*twiddleHandle)[i * 2]     = cos(M_PI * (forward ? -i : i) / (N));
         (*twiddleHandle)[i * 2 + 1] = sin(M_PI * (forward ? -i : i) / (N));
      }
      *lastNPointer = N;
   }
   return *twiddleHandle;
}


// In place complex input fft
// N is the number of complex values
// x should have 2*N values, interleaved with real and imaginary values
void rs_cfft(double* x, int N, int forward)
{
   // we use negative values with i when using the ifft
   int i;
   unsigned int j, twiddleIndex, width, blockStart, crawl, twiddleJump;
   double tempr, tempi;
   double scale = forward ? 1. / N : 1.;
   // reverse it
   rs_reverse(x, N);

   // this function guesses that we will be doing several fft calls of a given size
   // and caches the expensive sin values.
   double* twiddle;
   // for convenience use a pointer instead of a handle
   twiddle = fetch_twiddle_array(N, forward);

   // skip is the overall jump motion between groups,
   // width is the 
   for (width = 1, twiddleJump = N / 2; width < N; width *= 2, twiddleJump /= 2) {
      // our blocks are size of skip
      for (blockStart = 0; blockStart < N; blockStart += width * 2) {
         // walk through each of the subblocks with the skip
         // they go all the way to the end, skipping over more and more as
         // group size increases
         // this innermost loop jumps large bounds at the last bit

         // the way we jump over uneven blocks (e.g. second iteration in 8)
         // is that the blockStart determines the even/odd
         // and the inner loop just jumps skip
         for (crawl = 0; crawl < width; crawl++) {
            twiddleIndex = crawl * twiddleJump;
            i = blockStart + crawl;
            j = i + width;
            
            tempr = (x[j * 2]     * twiddle[twiddleIndex * 2] -
                     x[j * 2 + 1] * twiddle[twiddleIndex * 2 + 1]);
            tempi = (x[j * 2 + 1] * twiddle[twiddleIndex * 2] +
                     x[j * 2]     * twiddle[twiddleIndex * 2 + 1]);
            x[j * 2]     = x[i * 2]     - tempr;
            x[j * 2 + 1] = x[i * 2 + 1] - tempi;
            x[i * 2]     += tempr;
            x[i * 2 + 1] += tempi;
         }
      }
   }

   // not sure why but placing scaling before the FFT in the reverse 
   // function changes output, so we have to do it here.
   // perhaps we should use an inorder decimation in time, and reverse/scale after
   for (i = 0; i < N; i++) {
      x[i * 2]     *= scale;
      x[i * 2 + 1] *= scale;
   }
}

/* N should be set so that there are 2 * N real values in x */
void rs_rfft(double* x, int N, int forward)
{
   unsigned int i;
   double even_real, even_imag, odd_real, odd_imag;
   double *mirror;
   
   double* twiddle;

   if (forward) {
      rs_cfft(x, N, forward);
      // zero is a special case where we use the 0th complex value for both
      // the k and the mirror.
      // this makes the computation of twiddles a bit simpler too, since cos(0) = 1 
      // and sin(0) = 0

      even_real = x[0];
      // DC is just the sum of both 0 index vals
      x[0] = x[0] + x[1];
      // the N / 2 (nyquist uses the midpoint which also has special treatment.
      x[1] = even_real - x[1];
      // the midpoint - the real value is already set, but the other needs to be negated:
      x[N + 1] = -1. * x[N + 1];
   } else {
      even_real = x[0];
      // cancel out to get the original x[0]
      x[0] = (x[0] + x[1]) / 2;
      x[1] = (even_real - x[0]);
      x[N + 1] = -1 * x[N + 1];
   }

   // twiddle factors need to go from 0 to k * pi / N, but we can drop the pi/2 to pi second half -
   // we only need values from 1 to N / 2 because we use the symetrry around pi/2
   // fetch_twiddle_array gives us these values from 1 to N / 2.
   twiddle = fetch_twiddle_array_rfft(N, forward);

   for (i = 1; i < N / 2; i++) {
      // the formula for converting using symmetry is:
      // F(k) = .5 * ((Z(k) + Zconj(mirror(k))) - jW*(Z(k) -ZConj(mirror(k))))
      // where k is the complex bin index, and ZConj is the conjugate,
      // and mirror(k) is the symmetrical mirror index.
      // W is i * e^(-i*2*pi*k / N) where i is sqrt(-1) and N is the number of complex values

      // note that this is just a butterfly with the left side (Z(k) + Zconjf(mirror(k)))
      // being the even and the other side (- jW*(Z(k) -ZConj(mirror(k)))) being the odd.
      mirror = x + (N - i) * 2;

      // another way to look at the odd side is that after you factor the conjugate, 
      // it needs imaginary plus imaginary and real minus real,
      // and the even side is imaginary minus imaginary (because of the conjugate), and real plus real.
      even_real = 0.5 * (x[i * 2]     + mirror[0]);
      even_imag = 0.5 * (x[i * 2 + 1] - mirror[1]);
      odd_real  = (forward ? 0.5 : -0.5) * (x[i * 2]     - mirror[0]);
      odd_imag  = (forward ? 0.5 : -0.5) * (x[i * 2 + 1] + mirror[1]);

      // the e^(i*phi) = cos phi + i * sin phi splits the odd side up just like it does in a dft/fft 
      x[i * 2]     = even_real + (odd_real * twiddle[i * 2 + 1] + odd_imag * twiddle[i * 2]);
      x[i * 2 + 1] = even_imag - (odd_real * twiddle[i * 2]     - odd_imag * twiddle[i * 2 + 1]);
      // get the flip side as well.  All we have to do is reverse the twiddle sign for the real component
      // since cos is a negative mirror around pi / 2.
      // sin is a perfect mirror around pi / 2 so it may seem like we don't need to do anything.  However,
      // since the odd_real parts are subtracted, these need to be multiplied by -1 to flip it 
      // so that it is mirror - x[i * 2] instead.
      // so instead we can just add the odd part
      mirror[0]    = even_real - (odd_real * twiddle[i * 2 + 1]  + odd_imag * twiddle[i * 2]);
      // the flipping of even_imag applies just as it did with odd_real, so we subtract it.
      // the odd pairs don't need to be flipped because the imaginary component has
      // the subtraction and the cosine being multiplied, so the negatives cancel out
      mirror[1]    = -(odd_real * twiddle[i * 2]     - odd_imag * twiddle[i * 2 + 1]) - even_imag;
   }
   
   // the 0th index is special.
   // x[0] = 0.5;
   if (forward) {

   } else {
      rs_cfft(x, N, 0);
   }
}

/* uncomment to test 

#define NUMCOMPLEX 32

int main(int argc, char** argv)
{
   double x[NUMCOMPLEX * 2], y[NUMCOMPLEX * 2], z[NUMCOMPLEX * 2];
   double realx[NUMCOMPLEX], realy[NUMCOMPLEX], realz[NUMCOMPLEX];
   unsigned int i;

   // test reverse
   for (i = 0; i < NUMCOMPLEX; i++) {
      x[i * 2] = y[i * 2] = i;
      x[i * 2 + 1] = y[i * 2 + 1] = 0;
   }
   rs_reverse(x, NUMCOMPLEX);
   for (i = 0; i < NUMCOMPLEX; i++) {
      printf("x[%i]: %f, y[%i]: (%f, %f)\n", i, x[i * 2], i, y[i * 2], y[i * 2 + 1]);
   }

   printf("testing signal\n\n");

   // generate test sig
   for (i = 0; i < NUMCOMPLEX; i++) {
      x[i * 2] = y[i * 2] = realx[i] = realy[i] = (1 + // dc
                                                   cos(3 * i * 2 * M_PI / NUMCOMPLEX) + // 3 hz
                                                   sin(7 * i * 2 * M_PI / NUMCOMPLEX) + // near nyquist hz
                                                   sin(8 * i * 2 * M_PI / NUMCOMPLEX) + // special n/4 point
                                                   sin(12 * i * 2 * M_PI / NUMCOMPLEX) + // near nyquist hz
                                                   cos(i * M_PI) +
                                                   0); // nyquist
      x[i * 2 + 1] = y[i * 2 + 1] = 0;
   }
   rs_cfft(y, NUMCOMPLEX, 1);

   // copy to z, do inverse
   memcpy(z, y, sizeof(double) * NUMCOMPLEX * 2);
   rs_cfft(z, NUMCOMPLEX, 0);

   for (i = 0; i < NUMCOMPLEX; i++) {
      printf("x[%i]: %f, y[%i]: (%f, %f) z[%i]: (%f, %f)\n", i, x[i * 2], i, y[i * 2], y[i * 2 + 1], i, z[i * 2], z[i * 2 + 1]);
   }

   rs_rfft(realy, NUMCOMPLEX / 2, 1);
   // copy to z, do inverse
   memcpy(realz, realy, sizeof(double) * NUMCOMPLEX);
   rs_rfft(realz, NUMCOMPLEX / 2, 0);

   for (i = 0; i < NUMCOMPLEX; i++) {
      printf("x[%i]: %f, y[%i]: %f z[%i]: %f\n", i, realx[i], i, realy[i], i, realz[i]);
   }

   return 0;
}
*/

