#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/2
from lpp import *
from optparse import OptionParser
usage = '''usage: python2.7 %prog -a [ACE_file] -o [OUTPUT] -t [ Threshold  ]'''
parser = OptionParser(usage =usage ) 


parser.add_option("-f", "--Fgenesh", action="store", 
	dest="gene", 

	help="The Fgenesh File with mrna !!!!")

parser.add_option("-r", "--RAW", action="store", 
	dest="left", 
	default = 'left',
	help="The Fasta File input to Fgenesh !!!!")


(options, args) = parser.parse_args() 

Fegnesh = options.gene
fasta = options.left
FGENE_PREFIX = Fegnesh.split('.')[0]
FGENE = block_reading(  open( Fegnesh,'rU' )  ,tag = '//')
MRNA = open( FGENE_PREFIX+'.mrna' ,'w' )
PROTEIN = open( FGENE_PREFIX+'.pep' ,'w' )
DETAIL = open( FGENE_PREFIX+'.detail' ,'w' )
all_have = {}
OTHER_FASTA = open( FGENE_PREFIX+'.Error_fasta','w'  )
for each_b in FGENE:
	if ' no reliable predictions' in each_b:
		continue
	DETAIL.write( each_b+'\n//\n' )
	gene_name = re.search(  'Seq name\: (\S+) ',each_b ).group(1)
	all_have[ gene_name ] = ''
	seq_data = each_b.split( 'Predicted protein(s):\n' )[-1]
	all_gene_list = re.split( '\n(?=\>FGENESH\:\[mRNA\])',seq_data  )
	for each_data in all_gene_list:
		seq_list = re.split( '\n(?=\>)',each_data  )
		mrna  = seq_list[0]
		
		protein = seq_list[-1]
		number = re.search( '\s+(\d+)\s+', mrna).group(1)
		mrna_out = re.sub( '>FGENESH:\S*','>%s_%s'%( gene_name, number ),mrna+'\n'   )
		protein_out = re.sub( '>FGENESH:\S*','>%s_%s'%( gene_name, number ),protein+'\n'   )
		MRNA.write( mrna_out )
		PROTEIN.write(  protein_out  )
ALL_FASTA = fasta_check(  open( fasta,'rU'  )  )
END_FASTA = open( FGENE_PREFIX+'.F_fasta' ,'w' )
for t,s in ALL_FASTA:
	if t[1:-1] in all_have:
		END_FASTA.write(t+s)
	else:
		OTHER_FASTA.write( t+s  )