#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/20
from lpp import *
from optparse import OptionParser 
usage = '''usage: python2.7 %prog -o OUTPUT -a [Appendix]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output",
                  default = 'output', 
                  help="OUTPUT")
parser.add_option("-a", "--append", action="store", 
                  dest="append", 
                  default = 'contig',
                  help="The type you want")

(options, args) = parser.parse_args() 
output   = options.output
appendix = options.append
END = open(output+'.fasta','w')
LOG = open(output+'.list','w')
seq_all = Ddict()
for a,b,c in os.walk(os.getcwd()):
	for f in c:
		if f.endswith( '.%s'%( appendix ) ):
			
			RAW = fasta_check( open(a+'/'+f ,'rU'))
			
			for t,s in RAW:
				t= t[1:-1]
				s = re.sub( '\s+','',s  )
				if len(s) <=500:
					continue
				seq_all[s][f+'/'+t] = ''

i=0
LOG.write( 'ID\tFROM\n'  )
for each_seq in seq_all:
	i+=1
	END.write( '>%s\n'%( i )+each_seq+'\n' )
	LOG.write( '%s\t%s\n'%( i,'; '.join( seq_all[ each_seq ]  )   ) )
	