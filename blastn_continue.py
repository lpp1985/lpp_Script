#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/21
"""

from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU'))
all_seq = {}
for t,s in RAW:
	all_seq[t[1:-1]] = s.strip()+'\n'
BLAST_OUT = open(sys.argv[2],'rU')
data = BLAST_OUT.read()
all_have = re.findall("\<Iteration_query\-def\>([^<]+)",data)
already = {}
for key in all_have:
	already[key] = ""

while len(already)!=len(all_seq):
		
	QUERY = open("Query.pep",'w')
	for key,seq in all_seq.items():
		if key not in already:
			QUERY.write('>'+key+"\n"+seq)
	os.system(""" blastn -db /pub/Database/nt_2015_01_01/nt -query %s  -num_threads 63 -max_target_seqs 1 -evalue 1e-5  -outfmt  5    >> %s"""%(QUERY.name,sys.argv[2]  ))
	BLAST_OUT = open(sys.argv[2],'rU')
	data = BLAST_OUT.read()
	all_have = re.findall("\<Iteration_query\-def\>([^<]+)",data)
	already = {}
	for key in all_have:
		already[key] = ""	
