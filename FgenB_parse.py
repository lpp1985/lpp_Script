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

all_have = {}
ALL_FASTA = fasta_check(  open( fasta,'rU'  )  )

for t,s in ALL_FASTA:
	all_have[t[1:-1]] = re.sub( '\s+','', s  )
for each_b in FGENE:
	if ' no reliable predictions' in each_b:
		continue

	gene_name = re.search(  'Seq name\: (\S+) ',each_b ).group(1)
	seq_data = each_b.split( 'Predicted protein(s):\n' )[-1]
	all_gene_list = re.split( '\n(?=\>GENE)',seq_data  )
	for each_data in all_gene_list:
		data_list  = each_data.split( '\n' )
		title = data_list[0]

		pep = '\n'.join( data_list[1:] )+'\n'
		i_d = re.search( 'GENE\s+(\d+)',title  ).group(1)
		
		[start,stop] = [ int(x) for x in re.findall( '(\d+)\s+\-\s+(\d+)' ,title   )[0]]
		chain = re.search( 'chain\s+(\S)',title    ).group(1)
		mrna = all_have[gene_name  ][ start-1:stop  ]
		if chain =='-':
			mrna = complement( mrna )
		mrna = re.sub( '(\w{70})','\\1\n',mrna )
		name = '>%s_%s|%s:%s|chain_%s'%(  gene_name,i_d,start,stop,chain ) 
		MRNA.write( '%s\n%s\n'%(  name,mrna )   )
		PROTEIN.write( '%s\n%s\n'%( name,pep )  )

