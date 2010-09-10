
#ifndef ADBUFFER_LOADED
#define ADBUFFER_LOADED 1

extern void printtraffic() ;

extern void pushinteger4_(int x)  ;
extern void lookinteger4_(int *x) ;
extern void popinteger4_(int *x) ;

extern void pushreal4_(float x) ;
extern void lookreal4_(float *x) ;
extern void popreal4_(float *x) ;

extern void pushreal8_(double x) ;
extern void lookreal8_(double *x) ;
extern void popreal8_(double *x) ;

extern void printbuffertop() ;

#endif
