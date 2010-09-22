#ifndef DOT_HPP
#define DOT_HPP

inline int fid(int row, int col, int ldim){
    return row + col*ldim;
}

template <class Tdouble>
void dot(Tdouble *x, Tdouble *y, Tdouble *z, int N){
   for(int n = 0; n != N; ++n){ // iterate over all columns
       for(int m = 0; m != N; ++m){ // iterate over all rows
            z[fid(m,n,N)] = 0;
            for(int k = 0; k != N; ++k){
                z[fid(m,n,N)] += x[fid(m,k,N)]*y[fid(k,n,N)];
            }
        }
    }
}


#endif
