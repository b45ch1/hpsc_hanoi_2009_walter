import os
import sys
from pylab import *
import prettyplotting
import numpy
import re
import time


try:
	os.remove("runtimes.txt")
	os.remove("mem_consumption.txt")

except:
	pass


Ns = range(50,130,10)

for N in Ns:
    time.sleep(2)
    os.system("./matrix_ad_vs_adolc_vs_tapenade.exe %d" % N)

dat = numpy.loadtxt('runtimes.txt')

# Timings of the function evaluations
figure()
function_lapack = semilogy(dat[:,0], dat[:,1], '-k.')
function_my_c_implementation = semilogy(dat[:,0], dat[:,2], '-.k.')
# function_my_fortran_implementation = semilogy(dat[:,0], dat[:,3], 'k-.')
function_my_cpp_implementation = semilogy(dat[:,0], dat[:,4], '--k.')

adolc_taping   = semilogy(dat[:,0], dat[:,5], '-kv')
adolc_function =  semilogy(dat[:,0], dat[:,6], '--kv')


xlabel('matrix size $N$')
ylabel('runtime $t$ [ms]')

title(r'Function Evaluation')
legend((adolc_taping, adolc_function, function_lapack, function_my_c_implementation, function_my_cpp_implementation), ('ADOL-C (taping))', 'ADOL-C', 'LAPACK', 'Our C++ Impl.', 'Our C Impl.'),loc = 2)
grid()
savefig('function_matrix_ad_vs.eps')
savefig('function_matrix_ad_vs.png')

# Timing of the gradients
figure()
gradient_utps_adolc = semilogy(dat[:,0], dat[:,7], '-kv')
gradient_utpm_my_c_implementation = semilogy(dat[:,0], dat[:,8], '-.k.')
gradient_utpm_lapack = semilogy(dat[:,0], dat[:,9], '--k.')
gradient_utps_tapenade  = semilogy(dat[:,0], dat[:,10], '--kv')

xlabel('matrix size $N$')
ylabel('runtime $t$ [ms]')

title(r'Gradient Evaluation')
legend((gradient_utps_adolc, gradient_utps_tapenade, gradient_utpm_my_c_implementation,gradient_utpm_lapack), ('UTPS (ADOL-C)', 'UTPS (Tapenade)', 'UTPM (C)', 'UTPM (LAPACK)'),loc = 2)
grid()
savefig('gradient_matrix_ad_vs.eps')
savefig('gradient_matrix_ad_vs.png')

# Timing UTP arithmetic
figure()
utpm_taylorpoly = semilogy(dat[:,0], dat[:,-2], '-.k.')
utps_adolc = semilogy(dat[:,0], dat[:,-1], '-kv')
xlabel('matrix size $N$')
ylabel('runtime $t$ [ms]')

title(r'UTP Arithmetic')
legend((utpm_taylorpoly, utps_adolc), ('UTPM (TAYLORPOLY)', 'UTPS (ADOL-C)'),loc = 2)
grid()
savefig('utps_vs_utpm.eps')
savefig('utps_vs_utpm.png')


# #########################################
# #      MEMORY CONSUMPTION               #
# #########################################

# mem_list =[]


# p1 = re.compile('N=\d+')
# p = re.compile('\d+')
# f = open(r'mem_consumption.txt')
# N = 0
# for line in f.readlines():
# 	if len(p1.findall(line))>0:
# 		try:
# 			mem_list.append(tmplist)
# 		except:
# 			pass
# 		N = p.findall(line)[0]
# 		tmplist = [N]
# 		continue
# 	tmplist.append(p.findall(line)[0])
# 	#print p.findall(line)
# print asarray(mem_list)

# f.close()






show()