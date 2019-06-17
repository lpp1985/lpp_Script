#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2019/6/3
"""
from lpp import *
output = "./Input/"
os.mkdir(output)
for f in glob.glob("*.pep"):
	name = f.split(".")[0]
	PEP = open(output+name+".pep",'w')
	NUC = open(output+name+".nuc",'w')
	has_ = {}
	for t,s in fasta_check( open(f,'rU')):
		if "." in s:
			continue
		has_[t[1:].split()[0]]=""
		PEP.write(t+s)
	print(len(has_))
	for t,s in fasta_check( open(name+".nuc",'rU')):
		if t[1:].split()[0] in has_:
			del ( has_[t[1:].split()[0]])
			NUC.write(t+s)
	print( has_ )
