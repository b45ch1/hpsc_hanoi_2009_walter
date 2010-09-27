import os
import ctypes
import numpy
import time
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
lib.function_inv_cpp.argtypes = argtypes1
lib.function_inv_adolc.argtypes = argtypes1
lib.function_inv_atlas.argtypes = argtypes1


lib.gradient_trace_adolc.argtypes = argtypes1
lib.gradient_trace_c.argtypes = argtypes1
lib.gradient_trace_tapenade.argtypes = argtypes1
lib.gradient_trace_atlas.argtypes = argtypes1


lib.utp_dot_adolc.argtypes = argtypes2
lib.utp_dot_taylorpoly.argtypes = argtypes2
lib.utps_dot_taylorpoly.argtypes = argtypes2
runtimes = numpy.zeros(2, dtype=float)



N_list = [50,60,70,80,90,100,110,120,130,140,150]
# N_list = [40,50,]   


# # function evaluation
# num_methods = 5
# runtimes_matrix = numpy.zeros((num_methods, len(N_list)), dtype=float)
# for n,N in enumerate(N_list):
#     print 'N = ',N
#     runtimes[:] = 0
#     lib.function_inv_c(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[0,n] = runtimes[0]
    
#     runtimes[:] = 0
#     lib.function_inv_cpp(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[1,n] = runtimes[0]
    
#     runtimes[:] = 0
#     lib.function_inv_atlas(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[2,n] = runtimes[0]
    
#     runtimes[:] = 0
#     lib.function_inv_adolc(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[[3,4],n] = runtimes[:]
    
# pyplot.figure()
# pyplot.semilogy(N_list, runtimes_matrix[0,:], '-ko', markerfacecolor='None',label='Function eval (C)')
# pyplot.semilogy(N_list, runtimes_matrix[1,:], '--ko', markerfacecolor='None',label='Function eval (C++)')
# pyplot.semilogy(N_list, runtimes_matrix[2,:], '-kd', markerfacecolor='None',label='Function eval (ATLAS)')
# pyplot.semilogy(N_list, runtimes_matrix[3,:], '-.kv', markerfacecolor='None',label='ADOL-C tracing')
# pyplot.semilogy(N_list, runtimes_matrix[4,:], '--k^', markerfacecolor='None',label='ADOL-C eval')

# pyplot.xlabel('matrix size $N$')
# pyplot.ylabel('runtime $t$ [ms]')
# pyplot.title(r'Function Evaluation of inv(A)')
# pyplot.legend(loc=2)
# pyplot.grid()
# pyplot.savefig('function_evaluation_comparision.eps')

# N_list = [50,60,70,80,90,100,110,120,130,140,150]
# # gradient evaluation
# num_methods = 5
# runtimes_matrix = numpy.zeros((num_methods, len(N_list)), dtype=float)
# for n,N in enumerate(N_list):
#     print 'N = ',N
#     runtimes[:] = 0
#     lib.gradient_trace_c(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[0,n] = runtimes[0]
    
#     runtimes[:] = 0
#     lib.gradient_trace_atlas(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[1,n] = runtimes[0]
    
#     runtimes[:] = 0
#     lib.gradient_trace_adolc(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[[3,4],n] = runtimes[:]
    
#     runtimes[:] = 0
#     lib.gradient_trace_tapenade(N, runtimes.ctypes.data_as(c_double_ptr))
#     runtimes_matrix[2,n] = runtimes[0]



# pyplot.figure()
# pyplot.semilogy(N_list, runtimes_matrix[0,:], '-ko', markerfacecolor='None', label='UTPM eval (C)')
# pyplot.semilogy(N_list, runtimes_matrix[1,:], '-kd', markerfacecolor='None', label='UTPM eval (ATLAS)')
# pyplot.semilogy(N_list, runtimes_matrix[2,:], '-.k<', markerfacecolor='None', label='TAPENADE eval')
# pyplot.semilogy(N_list, runtimes_matrix[3,:], '-.kv', markerfacecolor='None',label='ADOL-C tracing')
# pyplot.semilogy(N_list, runtimes_matrix[4,:], '--k^', markerfacecolor='None',label='ADOL-C eval')
# pyplot.xlabel('matrix size $N$')
# pyplot.ylabel('runtime $t$ [ms]')
# pyplot.title(r'Gradient Evaluation trace(inv(A))')
# pyplot.legend(loc=2)
# pyplot.grid()
# pyplot.savefig('gradient_evaluation_comparision.eps')

# utp arithmetic
num_methods = 4
N_list = [50,60,70,80,90,100,110,120,130,140,150]
D_list = [1,2,3,10]
P_list = [1,5,10]
runtimes_matrix = numpy.zeros((num_methods, len(P_list), len(D_list), len(N_list)), dtype=float)

for p,P in enumerate(P_list):
    for d,D in enumerate(D_list):
        for n,N in enumerate(N_list):
            print 'P,D,N = ',P,D,N
            runtimes[:] = 0
            lib.utp_dot_taylorpoly(P, D, N, runtimes.ctypes.data_as(c_double_ptr))
            print runtimes
            runtimes_matrix[2,p,d,n] = runtimes[0]
            
            runtimes[:] = 0
            lib.utps_dot_taylorpoly(P, D, N, runtimes.ctypes.data_as(c_double_ptr))
            print runtimes
            runtimes_matrix[3,p,d,n] = runtimes[0]
            
            runtimes[:] = 0
            lib.utp_dot_adolc(P, D, N, runtimes.ctypes.data_as(c_double_ptr))
            print runtimes
            runtimes_matrix[0,p,d,n] = runtimes[0]
            runtimes_matrix[1,p,d,n] = runtimes[1]


for p,P in enumerate(P_list):
    for d,D in enumerate(D_list):
        pyplot.figure()
        pyplot.semilogy(N_list, runtimes_matrix[0,p,d,:], '-.kv', markerfacecolor='None', label='ADOL-C tracing')
        pyplot.semilogy(N_list, runtimes_matrix[1,p,d,:], '--k^', markerfacecolor='None', label='ADOL-C eval')
        pyplot.semilogy(N_list, runtimes_matrix[2,p,d,:], '-ko', markerfacecolor='None', label='TAYLORPOLY UTPM eval')
        pyplot.semilogy(N_list, runtimes_matrix[3,p,d,:], '-kd', markerfacecolor='None', label='TAYLORPOLY UTPS eval')
        pyplot.xlabel('matrix size $N$')
        pyplot.ylabel('runtime $t$ [ms]')
        pyplot.title(r'UTP, matrix-matrix product, P,D=%d,%d'%(P,D))
        pyplot.legend(loc=2)
        pyplot.grid()
        pyplot.savefig('utp_comparision_D%d_P%d.eps'%(D,P))

pyplot.show()

