#ifndef A6PulseShape_h
# define A6PulseShape_h

enum AtcPulseShapePar_t { ATCPAR_A0, ATCPAR_A, ATCPAR_T, ATCPAR_TAU, ATCPAR_P };

extern double A6PulseShape(double x,double par[5],double *der=0);

//for ROOT TF1 construction
extern double A6PulseShape(double *x,double *par);

#endif //A6PulseShape_h

