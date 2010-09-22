
import os
import ctypes
import numpy

lib = numpy.ctypeslib.load_library('libtest_cases', os.path.dirname(__file__))

c_double_ptr = ctypes.POINTER(ctypes.c_double)
c_int_ptr    = ctypes.POINTER(ctypes.c_int)
c_int        = ctypes.c_int
c_double     = ctypes.c_double

argtypes1 = [ctypes.c_int, c_double_ptr]
lib.function_inv_c.argtypes = argtypes1
lib.function_inv_adolc.argtypes = argtypes1
lib.gradient_trace_adolc.argtypes = argtypes1


runtimes = numpy.zeros(2, dtype=float)
for n in range(50,100,10):
    runtimes[:] = 0
    lib.function_inv_c(n, runtimes.ctypes.data_as(c_double_ptr))
    print runtimes
    runtimes[:] = 0
    lib.function_inv_adolc(n, runtimes.ctypes.data_as(c_double_ptr))
    print runtimes
    
    runtimes[:] = 0
    lib.gradient_trace_adolc(n, runtimes.ctypes.data_as(c_double_ptr))
    print runtimes
