#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

double* fetch_hamming(int n)
{
   static int lastN = -1;
   static double* hamming;
   if (lastN != n) {
      int i;
      float denom = 2. * M_PI / (n -1);
      lastN = n;

      if (hamming)
         free(hamming);
      hamming = (double*) malloc(sizeof(double) * n);
      for (i = 0; i < n; i++) {
         hamming[i] = .54 - .46 * cos(i * denom);
      }
   }
   return hamming;
}

// applies a hamming window to the input
void apply_hamming(double *input, int n)
{
   int i;
   double* hamming = fetch_hamming(n);
   for (i = 0; i < n; i++) {
      input[i] *= hamming[i];
   }
}
