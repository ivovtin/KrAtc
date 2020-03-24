#ifndef DFUNFT_H
# define DFUNFT_H

extern "C" void dfunft_(void sub(int*,double*,int*,double*,double*,double*,int*,int*),
	int *k, int *m, int *n, int *nx, int *nc, double *x, double *y, double *sy, double *a,
	double *al, double *au, int *mode, double *eps, int *maxit, int *iprt,
	int *mfr, int *iafr, double *phi, double *dphi, double *cov, double *std, double *w,
	int *nerr);

#endif
