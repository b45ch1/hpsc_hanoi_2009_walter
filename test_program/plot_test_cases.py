import os
import ctypes
import numpy
import matplotlib.pyplot as pyplot
import prettyplotting

lib = numpy.ctypeslib.load_library('libtest_cases', os.path.dirname(__file__))

c_double_ptr = ctypes.POINTER(ctypes.c_double)
c_int_ptr    = ctypes.POINTER(ctypes.c_int)
c_int        = ctypes.c_int
c_double     = ctypes.c_double

argtypes1 = [c_int, c_double_ptr]
argtypes2 = [c_int, c_int, c_int, c_double_ptr]

lib.function_inv_c.argtypes = argtypes1
lib.function_inv_adolc.argtypes = argtypes1
lib.gradient_trace_adolc.argtypes = argtypes1
lib.gradient_trace_c.argtypes = argtypes1

lib.utp_dot_adolc.argtypes = argtypes2
lib.utp_dot_taylorpoly.argtypes = argtypes2
lib.utps_dot_taylorpoly.argtypes = argtypes2





runtimes = numpy.zeros(2, dtype=float)

# # function evaluation
# for n in range(50,100,10):
#     print 'N = ',n
#     runtimes[:] = 0
#     lib.function_inv_c(n, runtimes.ctypes.data_as(c_double_ptr))
#     print runtimes
#     runtimes[:] = 0
#     lib.function_inv_adolc(n, runtimes.ctypes.data_as(c_double_ptr))
#     print runtimes

# gradient evaluation
# for n in range(50,130,10):
#     print 'N = ',n
#     runtimes[:] = 0
#     lib.gradient_trace_adolc(n, runtimes.ctypes.data_as(c_double_ptr))
#     print runtimes
    
#     runtimes[:] = 0
#     lib.gradient_trace_c(n, runtimes.ctypes.data_as(c_double_ptr))
#     print runtimes


# utp arithmetic

num_methods = 4
N_list = [50,80,100,120]
D_list = [10]
P_list = [1]

runtimes_matrix = numpy.zeros((num_methods, len(P_list), len(D_list), len(N_list)), dtype=float)

for p,P in enumerate(P_list):
    for d,D in enumerate(D_list):
        for n,N in enumerate(N_list):
            print 'P,D,N = ',P,D,N
            runtimes[:] = 0
            lib.utp_dot_adolc(P, D, N, runtimes.ctypes.data_as(c_double_ptr))
            runtimes_matrix[0,p,d,n] = runtimes[0]
            runtimes_matrix[1,p,d,n] = runtimes[1]
            runtimes[:] = 0
            lib.utp_dot_taylorpoly(P, D, N, runtimes.ctypes.data_as(c_double_ptr))
            runtimes_matrix[2,p,d,n] = runtimes[0]
            runtimes[:] = 0
            lib.utps_dot_taylorpoly(P, D, N, runtimes.ctypes.data_as(c_double_ptr))
            runtimes_matrix[3,p,d,n] = runtimes[0]

print runtimes_matrix
pyplot.figure()
pyplot.semilogy(N_list, runtimes_matrix[0,0,0,:], '-.kv', label='ADOL-C taping')
pyplot.semilogy(N_list, runtimes_matrix[1,0,0,:], '--k^', label='ADOL-C eval')
pyplot.semilogy(N_list, runtimes_matrix[2,0,0,:], '-k.', label='TAYLORPOLY UTPM eval')
pyplot.semilogy(N_list, runtimes_matrix[3,0,0,:], '-kd', label='TAYLORPOLY UTPS eval')

pyplot.xlabel('matrix size $N$')
pyplot.ylabel('runtime $t$ [ms]')
pyplot.title(r'UTP, matrix-matrix product, P,D=%d,%d'%(P_list[0],D_list[0]))
pyplot.legend(loc='best')

pyplot.show()

