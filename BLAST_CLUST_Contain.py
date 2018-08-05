#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/19
from lpp import *
from optparse import OptionParser 

usage = '''usage: python2.7 %prog -i input_path -t [The type you want]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  default = './', 
                  help="Input File")
parser.add_option("-c", "--Coverage", action="store", 
                  dest="coverage",
                  default = 0.85,
                  type = float, 
                  help="Input File")

parser.add_option("-o", "--end", action="store", 
                  dest="output", 
                  help="OUTPUT Data")

parser.add_option("-f", "--fasta", action="store", 
                  dest="fasta", 
                  help="Seq Data")
(options, args) = parser.parse_args()
if __name__ == '__main__':
	coverage = options.coverage
	INPUT = open(options.input, 'rU')
	END = open(options.output, 'w')
	all_has = {}
	all_length = {}
	for t, s in fasta_check( open(options.fasta)) :
		name = t.split()[0][1:]
		all_length[name] = len(s.replace("\n", ""))
	for line in INPUT:
		line_l = line.split()
		q_name = line_l[0]
		s_name = line_l[1]
		q_start= int(line_l[6]) 
		q_end =  int(line_l[7]) 
		q_start,q_end = sorted( [q_start,q_end]  )
		if q_name ==s_name:
			continue
		if q_name not in all_has:
			
			for q  in all_has:
				
				all_data = sorted(all_has[q], key = lambda x: x[0] )
				result_data = []
				for start, end in all_data:
					if not result_data:
						result_data.append([start, end])
					else:
						if start <= result_data[-1][-1] and end >= result_data[-1][-1] :
							result_data[-1][-1] = end
						elif start >result_data[-1][-1]:
							result_data.append( [start, end]  )
				aln_length = 0
				for start, end in result_data:
					aln_length += end - start + 1
				if float(aln_length) / all_length[ q] >= coverage:
					END.write("%s\t%d\t%d\n" % (q, all_length[ q] ,aln_length ))
			all_has = Ddict()
			all_has[q_name]= [[q_start, q_end]]
		else:
			all_has[q_name].append( [q_start, q_end])

	
		
	for q  in all_has:

		all_data = sorted(all_has[q], key = lambda x: x[0] )
		result_data = []
		for start, end in all_data:
			if not result_data:
				result_data.append([start, end])
			else:
				if start <= result_data[-1][-1] and end >= result_data[-1][-1] :
					result_data[-1][-1] = end
				elif start >result_data[-1][-1]:
					result_data.append( [start, end]  )
		aln_length = 0
		for start, end in result_data:
			aln_length += end - start + 1
		if float(aln_length) / all_length[ q] >= coverage:
			END.write("%s\t%d\t%d\n" % (q, all_length[ q] ,aln_length ))
	
