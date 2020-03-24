#ifndef CERNLIB_H
#define CERNLIB_H
// Этот файл нужен для правильного включения HBOOK в C-шные программы
typedef float  REAL;
typedef double REAL2;

/* MATHLIB */
extern float prob_( REAL*, int* );
extern void ranmar_( REAL*, int* );
extern float rndm_( REAL* );
//extern void ranlux_( REAL*, int* );
extern void rnormx_( REAL*, int*, int* );
extern void rnorml_( REAL*, int* );
extern void rnpssn_( REAL*, int*, int* );

/* HBOOK */
extern void hbook1_( int*, const char*, int*, REAL*, REAL*, REAL* , size_t );
extern void hbookn_( int*, const char*, int*, const char*, int*, const char*, size_t, size_t, size_t );
extern void hpak_( int*, REAL* );
extern void hunpak_( int*, REAL*, const char*, int*, size_t );
extern void hfill_( int*, REAL*, REAL*, REAL* );
extern void hfn_( int*, REAL* );
extern void hfithn_( int*, const char*, const char*, int*, REAL*, REAL*, REAL*, REAL*, REAL*, REAL*,size_t, size_t );
extern void hid1_( int*, int* );
extern void hid2_( int*, int* );
extern void hidall_( int*, int* );
extern void hdelet_( int* );
extern void hropen_( int*, const char*, const char*, const char*, int*, int*,
		    size_t, size_t, size_t );
extern void hrfile_( int*, const char*, const char*, size_t, size_t);
extern void hrout_( int*, int*, const char*, size_t );
extern void hrend_( const char*, size_t );
extern float rg32_( REAL* );
extern float hstati_( int*, int*, const char*, int*, size_t );
extern void igmeta_( int*, int*);
extern void igrng_( REAL*, REAL*);

/* HIGZ/HPLOT */
extern void iginit_( int* );
extern void igterm_( void );
extern void igend_( void );
extern void igsse_( int*, int* );
extern void igmeta_(int*, int*);
extern void hlimit_( int* );
extern void hplint_( int* );
extern void hplend_( void );
extern void igwkty_( int* );
extern void iopwk_( int*, int*, int* );
extern void iclwk_( int* );
extern void iclrwk_( int*, int* );
extern void iacwk_( int* );
extern void idawk_( int* );
extern void iuwk_( int*, int* );
extern void hplot_( int*, const char*, const char*, int*, size_t, size_t);
extern void hplzon_( int*, int*, int*, const char*, size_t );
extern void hplopt_( const char*, int*, size_t );
extern void hplset_( const char*, REAL*, size_t );
extern void igset_(const char*, REAL*, size_t );
extern void psfile_( int*, int*, const char*, size_t );

#endif /* CERNLIB_H */
