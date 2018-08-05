#!/usr/bin/python
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
parser.add_option("-i", "--input", action="store", 
                  dest="inputpath",
                  default = './', 
                  help="OUTPUT")
parser.add_option("-l", "--length", action="store", 
                  dest="length",
                  default = 400,
                   type= "int",
                  help="length")

parser.add_option("-n", "--name", action="store",
                  dest="name",
                  help="sequence name")

(options, args) = parser.parse_args() 

output   = options.output

length =  options.length
appendix = options.append

inputpath = options.inputpath
name = options.name
END = open(output+'.fasta','w')
LOG = open(output+'.list','w')
seq_all = Ddict()
for a,b,c in os.walk(os.path.abspath(inputpath)):
	for f in c:
		if f.endswith( '.%s'%( appendix ) ):
			try:
				RAW = fasta_check( open(a+'/'+f ,'rU'))
	
				for t,s in RAW:
					t= t[1:-1]
					s = re.sub( '\s+','',s  )
					if len(s)<=length:
						continue
					seq_all[s][f+': ' +t] = ''
			except:
				print( f )

i=0
LOG.write( 'ID\tFROM\n'  )
for each_seq in sorted(seq_all, key=lambda x: len(x) ):
	i+=1
	length = len( re.sub( '\s+','',each_seq    ) )
	END.write( '>%s_%s\n'%(name, i )+each_seq+'\n' )
	LOG.write( '%s_%s\t%s\n'%(name, i,'; '.join( seq_all[ each_seq ]  )   ) )
	
