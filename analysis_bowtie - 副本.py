#!/usr/bin/python
#Author=LPP
from lpp import *
import sys,re
#analysis_bowtie.py   [bowtie_end_file]  [  output_file   ]
nomis={}
onemis={}
twomis={}
no_mapping = Ddict()
one_mapping = Ddict()
two_mapping = Ddict()
INPUT = open( sys.argv[1] )
OUTPUT = open( sys.argv[2],'w' )
for line in INPUT:
	line_l = line[:-1].split('\t')
	mis = re.findall('(>)',line_l[-1])
	num = len( mis )
	try:
		if num>2:
			continue
		elif num ==2:
			two_mapping[ line_l[2] ][  line_l[0]  ]=''
			twomis[ line_l[0] ]=''
		elif num ==1:
			one_mapping[ line_l[2] ][  line_l[0]  ]=''
			onemis[ line_l[0] ]=''
		else:
			no_mapping[ line_l[2] ][  line_l[0]  ]=''
			nomis[ line_l[0] ]=''
	except:
		print(line)
		sys.exit()
for contig in no_mapping:
	OUTPUT.write(  contig +'\t'+'; '.join( no_mapping[ contig ] ) +'\t0\t\n' )
for contig in one_mapping:
	read = filter(lambda x: x not in nomis, one_mapping[ contig ])
	if read:
		OUTPUT.write(  contig +'\t'+'; '.join( read ) +'\t1\t\n' )
for contig in two_mapping:
	read = filter(lambda x: x not in nomis and x not in onemis, two_mapping[ contig ])
	if read:
		OUTPUT.write(  contig +'\t'+'; '.join( read ) +'\t2\t\n' )