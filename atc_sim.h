#ifndef _ATC_SIM_H
#define _ATC_SIM_H

#include "KcSys/fortran.h"

typedef struct {      //raw data from simulation
  int nHits;          //total number of hits
  float cou[500];     //number of counter
  float x1[500];
  float y1[500];
  float z1[500];
  float x2[500];
  float y2[500];
  float z2[500];
  float Ia_ch[500];
  float Ia_sc[500];
  float Ish_ch[500];
  float Ish_sc[500];
  float Itef_ch[500];
  float Itef_sc[500];
  float I_NAME[500];
  float I_NLEVEL[500];
} kscsim0_type;

extern kscsim0_type kscsim0_;

#ifdef __cplusplus
extern "C" {
#endif
    extern void atc_amp1_(int* a1, float* a2,float* a3);
    extern void atc_amp_coef_(int* x1, float* x2,float* x3, float* x4,float* x5, float* x6,float* x7, float* x8,int* x9,int* x10);
    extern void atc_hist_();
#ifdef __cplusplus
}
#endif


#endif


