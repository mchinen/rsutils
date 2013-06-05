/*
 *  rsutil.h
 *  Filter
 *
 *  Created by mchinen on 11/23/10.
 *  Copyright 2010 Roughsoft. All rights reserved.
 *
 */

#ifndef __RSUTIL_H__
#define __RSUTIL_H__

#include "stdint.h"

#ifdef __cplusplus
extern "C" {
#endif

#define TWOPI 6.283185
#define PI    3.141592
//macros

//when using SIMP_MAX be sure not to have any functions that have side effects
//as parameters or it will blow up as it gets evaluated twice.
//use funcmax() for that case.
#define SIMP_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SIMP_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SIMP_CLAMP(val,min,max) ( ((val) < (min)) ? (min) : ( ((val) > (max)) ? (max) : (val) ) ) 


double eval_min(double a, double b);
double eval_max(double a, double b);
double eval_clamp(double val, double min, double max);
//rolls a n-sided die that starts at zero
int randint(int n);
//float under specified value
float randfloat(float under);


double amptodb(double);
double dbtoamp(double);

double sumsquared(double pair[2]);
double rootsumsquared(double pair[2]);

double round(double r);

float randfloatexp2(float under, int divs);
float exp2cursor(float cursor);
int coinflip();

// super basic checksum from two uints.
// you can reuse this to make more complicated ones.
uint32_t ezchecksum(uint32_t a, uint32_t b);

#ifdef __cplusplus
};
#endif

#endif
