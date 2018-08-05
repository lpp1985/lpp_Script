#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/3
"""

from lpp import *

from optparse import OptionParser
usage = '''usage: python2.7 %prog [options] '''
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Input", action="store",
                             dest="Input",
)
parser.add_option("-u", "--Unigene", action="store",
                      dest="Unigene",
                      help="Unigene(Longest isoform)")
parser.add_option("-l", "--List", action="store",
                      dest="List",
                      help="List of Unigene")

(options, args) = parser.parse_args()
Input = options.Input
Unigene = open(options.Unigene,'w')
List = open(options.List,'w')

RAW = fasta_check(open(Input,'rU'))
data_hash = Ddict()
List.write("Gene\tIsoforms\tIsoform.No\n")
for t,s in RAW:
	s1 = re.sub("\s+","",s)
	name = t[1:].split()[0]
	gene_name = name.rsplit("_",1)[0]
	data_hash[gene_name][len(s1)][name] = s
for e_geneName  in data_hash:
	all_seq = []
	iso_name = sorted(data_hash[e_geneName])[-1]
	for name,seq in data_hash[e_geneName][iso_name].items():
		Unigene.write('>'+name+'\n'+seq)
		break
	for length in data_hash[e_geneName]:
		for name in data_hash[e_geneName][length]:
			all_seq.append(name)
	List.write(e_geneName+'\t'+'; '.join(all_seq)+'\t%s\n'%(len(all_seq) ) )