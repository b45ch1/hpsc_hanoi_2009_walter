static char adSid[]="$Id: adBuffer.c $";

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "adStack.h"

/************ MEASUREMENT OF PUSH/POP TRAFFIC *************/

static int mmftraffic = 0 ;
static int mmftrafficM = 0 ;

void addftraffic(int n) {
  mmftraffic = mmftraffic+n ;
  while (mmftraffic >= 1000000) {
    mmftraffic = mmftraffic-1000000 ;
    ++mmftrafficM ;
  }
  while (mmftraffic < 0) {
    mmftraffic = mmftraffic+1000000 ;
    --mmftraffic ;
  }
}

void printtraffic() {
  printctraffic_() ;
  printf(" F Traffic: ") ;
  printbigbytes(mmftrafficM, 1000000, mmftraffic) ;
  printf(" bytes\n") ;
}

/************************** integer*4 ************************/
static int adi4buf[512] ;
static int adi4ibuf = 0 ;
static int adi4lbuf[512] ;
static int adi4ilbuf = -1 ;
static int adi4inlbuf = 0 ;

void pushinteger4_(int x) {
  addftraffic(4) ;
  if (adi4ilbuf != -1) {
    adi4ilbuf = -1 ;
    adi4inlbuf = 0 ;
  }
  if (adi4ibuf >= 511) {
    adi4buf[511] = x ;
    pushNarray((char*)adi4buf, 512*4) ;
    addftraffic(-512*4) ;
    adi4ibuf = 0 ;
  } else {
    adi4buf[adi4ibuf] = x ;
    ++adi4ibuf ;
  }
}

void lookinteger4_(int *x) {
  if (adi4ilbuf == -1) {
    adi4ilbuf = adi4ibuf ;
    resetadlookstack_() ;
  }
  if (adi4ilbuf <= 0) {
    lookNarray((char*)adi4lbuf, 512*4) ;
    adi4inlbuf = 1 ;
    adi4ilbuf = 511 ;
    *x = adi4lbuf[511] ;
  } else {
    --adi4ilbuf ;
    if (adi4inlbuf)
      *x = adi4lbuf[adi4ilbuf] ;
    else
      *x = adi4buf[adi4ilbuf] ;
  }
}

void popinteger4_(int *x) {
  if (adi4ilbuf != -1) {
    adi4ilbuf = -1 ;
    adi4inlbuf = 0 ;
  }
  if (adi4ibuf <= 0) {
    popNarray((char*)adi4buf, 512*4) ;
    adi4ibuf = 511 ;
    *x = adi4buf[511] ;
  } else {
    --adi4ibuf ;
    *x = adi4buf[adi4ibuf] ;
  }
}

/************************** real*4 ************************/
static float adr4buf[512] ;
static int adr4ibuf = 0 ;
static float adr4lbuf[512] ;
static int adr4ilbuf = -1 ;
static int adr4inlbuf = 0 ;

void pushreal4_(float x) {
  addftraffic(4) ;
  if (adr4ilbuf != -1) {
    adr4ilbuf = -1 ;
    adr4inlbuf = 0 ;
  }
  if (adr4ibuf >= 511) {
    adr4buf[511] = x ;
    pushNarray((char*)adr4buf, 512*4) ;
    addftraffic(-512*4) ;
    adr4ibuf = 0 ;
  } else {
    adr4buf[adr4ibuf] = x ;
    ++adr4ibuf ;
  }
}

void lookreal4_(float *x) {
  if (adr4ilbuf == -1) {
    adr4ilbuf = adr4ibuf ;
    resetadlookstack_() ;
  }
  if (adr4ilbuf <= 0) {
    lookNarray((char*)adr4lbuf, 512*4) ;
    adr4inlbuf = 1 ;
    adr4ilbuf = 511 ;
    *x = adr4lbuf[511] ;
  } else {
    --adr4ilbuf ;
    if (adr4inlbuf)
      *x = adr4lbuf[adr4ilbuf] ;
    else
      *x = adr4buf[adr4ilbuf] ;
  }
}

void popreal4_(float *x) {
  if (adr4ilbuf != -1) {
    adr4ilbuf = -1 ;
    adr4inlbuf = 0 ;
  }
  if (adr4ibuf <= 0) {
    popNarray((char*)adr4buf, 512*4) ;
    adr4ibuf = 511 ;
    *x = adr4buf[511] ;
  } else {
    --adr4ibuf ;
    *x = adr4buf[adr4ibuf] ;
  }
}

/************************** real*8 ************************/
static double adr8buf[512] ;
static int adr8ibuf = 0 ;
static double adr8lbuf[512] ;
static int adr8ilbuf = -1 ;
static int adr8inlbuf = 0 ;

void pushreal8_(double x) {
  addftraffic(8) ;
  if (adr8ilbuf != -1) {
    adr8ilbuf = -1 ;
    adr8inlbuf = 0 ;
  }
  if (adr8ibuf >= 511) {
    adr8buf[511] = x ;
    pushNarray((char*)adr8buf, 512*8) ;
    addftraffic(-4096) ;
    adr8ibuf = 0 ;
  } else {
    adr8buf[adr8ibuf] = x ;
    ++adr8ibuf ;
  }
}

void lookreal8_(double *x) {
  if (adr8ilbuf == -1) {
    adr8ilbuf = adr8ibuf ;
    resetadlookstack_() ;
  }
  if (adr8ilbuf <= 0) {
    lookNarray((char*)adr8lbuf, 512*8) ;
    adr8inlbuf = 1 ;
    adr8ilbuf = 511 ;
    *x = adr8lbuf[511] ;
  } else {
    --adr8ilbuf ;
    if (adr8inlbuf)
      *x = adr8lbuf[adr8ilbuf] ;
    else
      *x = adr8buf[adr8ilbuf] ;
  }
}

void popreal8_(double *x) {
  if (adr8ilbuf != -1) {
    adr8ilbuf = -1 ;
    adr8inlbuf = 0 ;
  }
  if (adr8ibuf <= 0) {
    popNarray((char*)adr8buf, 512*8) ;
    adr8ibuf = 511 ;
    *x = adr8buf[511] ;
  } else {
    --adr8ibuf ;
    *x = adr8buf[adr8ibuf] ;
  }
}

/********* PRINTING THE SIZE OF STACKS AND BUFFERS ********/

void printbuffertop() {
  int size = 0 ;
  size += adi4ibuf*4 ;
  size += adr4ibuf*4 ;
  size += adr8ibuf*8 ;
  printf("Buffer size:%i bytes i.e. %i Kbytes\n",
         size, size/1024.0) ;
}

/**********************************************************
 *        HOW TO CREATE PUSH* LOOK* POP* SUBROUTINES
 *              YET FOR OTHER DATA TYPES
 * Duplicate and uncomment the commented code below.
 * In the duplicated code, replace:
 *   ctct -> C type name (e.g. float double, int...)
 *   TTTT -> BASIC TAPENADE TYPE NAME
 *     (in character, boolean, integer, real, complex, pointer,...)
 *   z7   -> LETTER-SIZE FOR TYPE
 *     (in s,         b,       i,       r,    c,       p,      ...)
 *   7    -> TYPE SIZE IN BYTES
 * Don't forget to insert the corresponding lines in
 * procedure printbuffertop(), otherwise the contribution of
 * this new type to buffer occupation will not be seen.
 * (not very important anyway...)
 **********************************************************/

/************************** TTTT*7 ************************/
/*
static ctct adz7buf[512] ;
static int adz7ibuf = 0 ;
static ctct adz7lbuf[512] ;
static int adz7ilbuf = -1 ;
static int adz7inlbuf = 0 ;

void pushTTTT7_(ctct x) {
  addftraffic(7) ;
  if (adz7ilbuf != -1) {
    adz7ilbuf = -1 ;
    adz7inlbuf = 0 ;
  }
  if (adz7ibuf >= 511) {
    adz7buf[511] = x ;
    pushNarray((char*)adz7buf, 512*7) ;
    addftraffic(-512*7) ;
    adz7ibuf = 0 ;
  } else {
    adz7buf[adz7ibuf] = x ;
    ++adz7ibuf ;
  }
}

void lookTTTT7_(ctct *x) {
  if (adz7ilbuf == -1) {
    adz7ilbuf = adz7ibuf ;
    resetadlookstack_() ;
  }
  if (adz7ilbuf <= 0) {
    lookNarray((char*)adz7lbuf, 512*7) ;
    adz7inlbuf = 1 ;
    adz7ilbuf = 511 ;
    *x = adz7lbuf[511] ;
  } else {
    --adz7ilbuf ;
    if (adz7inlbuf)
      *x = adz7lbuf[adz7ilbuf] ;
    else
      *x = adz7buf[adz7ilbuf] ;
  }
}

void popTTTT7_(ctct *x) {
  if (adz7ilbuf != -1) {
    adz7ilbuf = -1 ;
    adz7inlbuf = 0 ;
  }
  if (adz7ibuf <= 0) {
    popNarray((char*)adz7buf, 512*7) ;
    adz7ibuf = 511 ;
    *x = adz7buf[511] ;
  } else {
    --adz7ibuf ;
    *x = adz7buf[adz7ibuf] ;
  }
}
*/

/**********************************************************/
