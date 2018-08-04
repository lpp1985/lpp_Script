#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
# Company: Chinese Human Genomic Center at Beijing
from lpp import *
from optparse import OptionParser
usage='''usage: python %prog [options] '''
parser = OptionParser(usage =usage )
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input_fasta",
                  type='string',  
                  help="input fasta (scaffold file)")
parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="outputDIR",
                  type='string',  
                  help="OUTPUT DIR")

(options, args) = parser.parse_args() 
if __name__=='__main__':
	outpurt_dir = os.path.abspath( options.outputDIR  )+'/'
	RAW = fasta_check( open( options.input_fasta,'rU'  )   )
	CONTIG = open( outpurt_dir+'contigs.fa.original','w'      )
	ORG = open( outpurt_dir+'read.placed.original','w'   )
	for a,b in RAW:
		b = re.sub( '\s+','',b  )
		b_block = re.split( 'N+',b )
		print( len(  b_block) )
		i=0
		title = a[1:-1]
		for each_b in b_block:
			if len(  each_b  )<100:
				continue
			i+=1
			name = title+'_fragment%s'%( i )
			ORG.write( name+'\t'+ title+'\n' )
			CONTIG.write( '>'+name+'\n'+each_b+'\n'  )