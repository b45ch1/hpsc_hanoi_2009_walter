#!/usr/bin/env python
import os

#try to compile first
os.system('g++ -c dense_inverse.c givens.c')
# os.system('clear')

# differentiate
os.system('tapenade -reverse -head inv -vars "a" -outvars "qt" -inputlanguage C -outputlanguage C -html dense_inverse.c')
os.system('tapenade -reverse -head givens -vars "a b" -outvars "c s" -inputlanguage C -outputlanguage C -html givens.c')

