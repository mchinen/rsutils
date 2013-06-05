/*
 *  rsutil.c
 *  Filter
 *
 *  Created by mchinen on 11/23/10.
 *  Copyright 2010 Roughsoft. All rights reserved.
 *
 */

#include "math.h"
#include "stdlib.h"

#include "rs_util.h"

// if you want to use a custom mask for our crappy checksum, you can do this in
// your compiler flags with a -D
#ifndef EZCHECKSUMMASK
#define EZCHECKSUMMASK 0xF0F1
#endif

int randint(int n)
{
    return rand() % n;
}

float randfloat(float under)
{
    return (rand() / (RAND_MAX + 1.0)) * under;
}

double eval_min(double a, double b)
{
    return SIMP_MIN(a, b);
}
double eval_max(double a, double b)
{
    return SIMP_MAX(a, b);
}

double eval_clamp(double val, double min, double max)
{
   return SIMP_CLAMP(val, min, max);
}

double amptodb(double amp)
{
   return 20 * log10((double)amp);
}

double dbtoamp(double db)
{
   return powf(10,  db/20.0);
}

inline double sumsquared(double pair[2])
{
   return  pair[0] * pair[0] + pair[1] * pair[1];
}

inline double rootsumsquared(double pair[2])
{
   return  sqrt(pair[0] * pair[0] + pair[1] * pair[1]);
}

inline double round(double r)
{
   return r > 0.0 ? floor(r + 0.5) : ceil(r - 0.5);
}


//expontential in div by 2 - the last 50% has the same chance as the 25-50%,12.5-25%bracket, and so on.  this will make
//lower values more likely, but it can't go on for infinity.  we cut it once we reach a 1% bracket.

//this next number essentially should be the number of octaves you want.
#define kRandFloatExp2Limit 9 
float randfloatexp2(float under, int divs)
{
   int part = randint(kRandFloatExp2Limit);
   
   if(part>0) 
      return (1.0/(1 << part) + randfloat(1.0/(1<<part)))*under;
   else
      {
         //special case for the remainder.
         return randfloat(1.0/(1<<((divs<1?kRandFloatExp2Limit:divs)-1)));
      }
}

//cursor from 0 to 1, return a logarithmic cursor s.t. the 50%-100% bracket has the same cursor width as the 25-50%etc.  
//we use randfloatexp2limit to define the smallest unit.
float exp2cursor(float cursor)
{
   if(cursor<2.0/((float)kRandFloatExp2Limit))
      return cursor*kRandFloatExp2Limit/2.0 * powf(2.0,(2.0/((float)kRandFloatExp2Limit) - 1.0)* kRandFloatExp2Limit);
   else
      return powf(2.0,(cursor - 1.0)* kRandFloatExp2Limit);
}

int coinflip()
{
   return randint(2);
}

uint32_t ezchecksum(uint32_t a, uint32_t b)
{
   return a + ((~a) & EZCHECKSUMMASK) + (b & 5) + (a ^ b);
}
