#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/6/30
"""

from lpp import *
RAW = open(sys.argv[1],'rU')
ref_has = {}
seq_has = {}
for line in RAW:
	line_l = line.split("\t")
	for key in line_l[2].strip().split()[:-1]:
		ref_has[key] = ''
	for key in line_l[-1].strip().split()[:-1]:
		seq_has[key] = ''
DEL = open("Delete.fasta",'w')
DEL_LIST = open("Deletion.list",'w')
DEL_LIST.write("Gene\tAnnotation")
for t,s in fasta_check(open(sys.argv[2],'rU')):
	name = t[1:].split()[0]
	if name not in ref_has:
		DEL.write(t+s)
		DEL_LIST.write("\t".join(t.strip()[1:].split())+'\n')
		
		
UNIQUE = open("Unique.fasta",'w')
UNI_LIST = open("Unique.list",'w')
UNI_LIST.write("Gene\tAnnotation")
for t,s in fasta_check(open(sys.argv[3],'rU')):
	name = t[1:].split()[0]
	if name not in seq_has:
		UNIQUE.write(t+s)
		UNI_LIST.write("\t".join(t.strip()[1:].split())+'\n')	
