#include <math.h>
int givens(double *a, double *b, double *c, double *s){
    /*
    Algorithm that computes c and s from a and b s.t.
    
      [ c   s  ] [a]     [r]
      [ -s   c ] [b]  =  [0]
    
    is satisfied.
    */
    double r;
    if(*b == 0.){
        *c = 1; *s = 0;
    }
    else{
        if(fabs(*b) > fabs(*a)){
             r = - (*a)/(*b);
            *s = 1./sqrt(1 + r*r);
            *c = *s * r;
        }
        else{
             r = - (*b)/(*a);
            *c = 1./sqrt(1 + r*r);
            *s = *c * r;
        }
    }
}
