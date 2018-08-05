#!/usr/bin/env python
#Author=LPP
from lpp import *
import subprocess, string, shutil,multiprocessing
END= open( 'ALL_gap_static.fasta' ,'w')

RAW = fasta_check(open( sys.argv[1] ))
TEMP = open(  'Template.fasta','w' )
STA = open(  'ALL_gap.location','w' )
for t,s in RAW:
	title = re.search( '(\d+)',t ).group(1)
	s1 = re.sub('\s+','',s)
	i=1
	TEMP.write( t+s1+'\n' )
	for each_Gap in re.finditer('(N+)',s1,re.I):
		[ gap_start, gap_stop ] = [ int(x) for x in each_Gap.span() ]
		STA.write( '>s'+title+'-%s'%(i)+'\t'+'%s\t%s\t%s'%( gap_start, gap_stop,gap_stop-gap_start )  +'\n' )
		[ gap_start, gap_stop ] = [  gap_start-150,gap_stop + 150    ]
		END.write( '>s'+title+'-%s'%(i)+'\n'+s1[ gap_start: gap_stop ]+'\n' )
		i=i+1
