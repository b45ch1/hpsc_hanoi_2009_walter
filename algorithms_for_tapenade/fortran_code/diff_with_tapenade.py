#!/usr/bin/env python
import os

#try to compile first
os.system('g++ -c dense_inverse.c')
os.system('clear')

# differentiate
os.system('tapenade -reverse -head inv -vars "a" -outvars "qt" -inputlanguage fortran -outputlanguage fortran -html dense_inverse.f dense_qr_decomposition_by_givens_rotations.f')
 
