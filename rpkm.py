#!/usr/bin/env python
#coding:utf-8
"""
Author:   --<>
Purpose: 
Created: 2014/11/3
"""

import numpy as np
#import matplotlib.pyplot as plt
import re
from lpp import *
import glob
from  optparse import OptionParser
usage='''usage: python %prog [options]
It can automaticaly analysis RPKM and plot it'''

parser = OptionParser(usage =usage )
#parser.add_option("-i", "--Input", action="store",
                  #dest="input",
                  #type='string',
                  #help="input fastq file")
parser.add_option("-f", "--fasta_file", action="store",
                  dest="fasta",
                  type='string',

                  help="unigene fasta ")
parser.add_option("-m", "--Matrix_file", action="store",
                  dest="matrix",
                  type='string',

                  help="reads mapping number")


parser.add_option("-o", "--OUTPUT", action="store",
                  dest="output",
                  type='string',

                  help="RPKM output")
parser.add_option("-g", "--Graph", action="store",
                  dest="graph",
                  type='string',

                  help="RPKM graph  output")


(options, args) = parser.parse_args()
#input_path = options.input

#all_fastq_file = glob.glob(input_path+"/*.pair1")
#data_file_group = {}
#for each_file in all_fastq_file:
	#name = os.path.split(each_file)[-1].split('.')[0].split("_")[0]
	#data_file_group[name] = each_file
   
reads_count = {}
MATRIX = open(options.matrix,'rU')
title = MATRIX.next()
title_l = title.strip().split("\t")
for key in title_l[1:]:
	reads_count[key] = 0
for line in MATRIX:
	line_l = line.split("\t")
	for i in xrange(1, len(line_l)):
		reads_count[title_l[i]] += int( line_l[i])



#for name, data  in data_file_group.items():
	#line_num = int(os.popen("wc -l %s"%(data)).read().split()[0])/2
	#reads_count[name] = line_num

unigene_length = {}
RAW = fasta_check(open(options.fasta))
for t,s in RAW:
	unigene_length[t[1:].split()[0]] = len(re.sub("\s+", '', s))
END = open(options.output,'w')
MATRIX = open(options.matrix,'rU')



title = MATRIX.next()
title_l = title.strip().split("\t")[1:]
data_hash = {}
all_data = []
for each_sample in title_l:
	data_hash[each_sample] = []
END.write(title)
def rpkm_trans(number):
	number = float(number)

for line in MATRIX:
	line_l = line.strip().split("\t")
	gene_name = line_l[0]
	express_l = line_l[1:]
	rpkm_cache = []
	for i in xrange(0,len(express_l) ):
		rpkm = 10**9*float(express_l[i])/reads_count[title_l[i]]/unigene_length[gene_name]
		all_data.append(rpkm)
		data_hash[ title_l[i]    ].append(float(express_l[i]) )
		rpkm_cache.append("%s"%(rpkm))
	END.write(gene_name+"\t"+"\t".join(rpkm_cache)+'\n')
		
#data = [ data_hash[i]  for i in sorted(data_hash)]
#fig, ax1 = plt.subplots(figsize=(10,6),dpi=100)
#fig.canvas.set_window_title('Expression Level')
#plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

#bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
#plt.setp(bp['boxes'], color='black')
#plt.setp(bp['whiskers'], color='black')
#plt.setp(bp['fliers'], color='red', marker='+')

## Add a horizontal grid to the plot, but make it very light in color
## so we can use it for reading data values but not be distracting
#ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)
#ax1.set_yscale("log")
# Hide these grid behind plot objects
#ax1.set_axisbelow(True)
#ax1.set_title('Different Sample\'s Expression Value')
#ax1.set_xlabel('Distribution')
#ax1.set_ylabel('RPKM')
#xtickNames = plt.setp(ax1, xticklabels=[ "sample %s"%(i) for i in  sorted(data_hash) ])
#ax1.set_ylim(  (0.000001,max(all_data))  )  

#plt.savefig("%s_RPKM.png"%(options.graph))
#plt.show()
