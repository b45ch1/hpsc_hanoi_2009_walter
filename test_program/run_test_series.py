import os
import sys
from pylab import *
import numpy
import re

try:
	os.remove("runtimes.txt")
	os.remove("mem_consumption.txt")

except:
	pass


Ns = range(30,160,10)

for N in Ns:
	os.system("matrix_ad_vs_adolc_vs_tapenade.exe %d" % N)



dat = numpy.load('runtimes.txt')


# Timings of the function evaluations
figure()
adolc_taping   = semilogy(dat[:,0], dat[:,1], 'k.')
adolc_function =  semilogy(dat[:,0], dat[:,7], 'k--')
atlas_implementation = semilogy(dat[:,0], dat[:,2], 'k:')
my_c_implementation = semilogy(dat[:,0], dat[:,3], 'k-')
my_F77_implementation = semilogy(dat[:,0], dat[:,4], 'k-.')

xlabel('matrix size $N$')
ylabel('runtime $t$ [ms]')

title(r'Function Evaluation')
legend((adolc_taping, adolc_function, atlas_implementation, my_c_implementation, my_F77_implementation), ('ADOL-C (taping))', 'ADOL-C', 'Atlas', 'Our C++ Impl.', 'Our F77 Impl.'),loc = 2)
grid()
savefig('function_matrix_ad_vs.eps')
savefig('function_matrix_ad_vs.png')

# Timing of the gradients
figure()
adolc_plot     = semilogy(dat[:,0], dat[:,8], 'k--')
matrix_ad_plot = semilogy(dat[:,0], dat[:,9], 'k-')
tapenade_plot  = semilogy(dat[:,0], dat[:,10], 'k:')



xlabel('matrix size $N$')
ylabel('runtime $t$ [ms]')

title(r'UTPS vs UTPM')
legend((adolc_plot,matrix_ad_plot,tapenade_plot), ('UTPS (ADOL-C)', 'UTPM', 'UTPS (Tapenade)'),loc = 2)
grid()
savefig('gradient_matrix_ad_vs.eps')
savefig('gradient_matrix_ad_vs.png')

#########################################
#      MEMORY CONSUMPTION               #
#########################################

mem_list =[]


p1 = re.compile('N=\d+')
p = re.compile('\d+')
f = open(r'mem_consumption.txt')
N = 0
for line in f.readlines():
	if len(p1.findall(line))>0:
		try:
			mem_list.append(tmplist)
		except:
			pass
		N = p.findall(line)[0]
		tmplist = [N]
		continue
	tmplist.append(p.findall(line)[0])
	#print p.findall(line)
print asarray(mem_list)

f.close()






#show()