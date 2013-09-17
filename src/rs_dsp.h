/****************************************************************
 
    rs_dsp.h

    DSP utility calls

    Michael Chinen
    http://roughsoft.com

    Copyright 2012 Roughsoft LLC
    License: MIT

****************************************************************/

#ifndef __RS_DSP_H__
#define __RS_DSP_H__


#ifdef __cplusplus
extern "C" {
#endif


double* fetch_hamming(int n);
void apply_hamming(double *input, int n);

#ifdef __cplusplus
};
#endif

#endif
