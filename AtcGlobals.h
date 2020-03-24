#ifndef AtcGlobals_h
# define AtcGlobals_h 1

/*Number of all ATC counters*/
# ifndef ATC_NCNT
#  define ATC_NCNT 160
# endif

/*Number of installed ATC counters*/
# ifdef NATC
#  undef NATC
# endif
static const int NATC=160;

/*Maximum number of stored DC tracks*/
# ifndef ATC_NTRK
#  define ATC_NTRK 16
# endif

/*Maximum number of stored counters bound to track*/
# ifndef ATC_NCROS
#  define ATC_NCROS 16
# endif

#endif
