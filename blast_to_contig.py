#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
# Company: Chinese Human Genomic Center at Beijing
from lpp import *
from optparse import OptionParser
#from quality_stats_singleEnd import fastq_quality_class
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
usage='''usage: python %prog [options] '''
parser = OptionParser(usage =usage )
parser.add_option("-b", "--blast", action="store", 
                  dest="blast",
                  type='string',  
                  help="Blast_Result")
parser.add_option("-s", "--scaffold", action="store", 
                  dest="scaffold",
                  type='string',  
                  help="Scafoold file(fasta format)")
parser.add_option("-r", "--reference", action="store", 
                  dest="reference",
                  type='string',  
                  help="cosmid 2 reference file")
parser.add_option("-o", "--output", action="store", 
                  dest="output",
                  type='string',  
                  help="output contig file")

(options, args) = parser.parse_args() 
PAIR= open( 'test.mates','w'  )
BLAST = open(  options.blast , 'rU' )
REFBLAST = open( options.reference ,'rU'  )
SCAFFOLD = fasta_check( open(  options.scaffold , 'rU' ) )
all_seq = {}
OUTPUT = open(  options.output ,'w' )
OTHER = open( 'other.scaff','w'  )
all_cosmid = Ddict()
for t,s in SCAFFOLD:
	s = re.sub( '\s+','',s )
	title = t.split()[0][1:]
	all_seq[ title ] = s

if __name__=='__main__':
	
	all_blast_detail = {}  
	no_unique = {}
	BLAST.next()
	for line in BLAST:
		line_l = line[:-1].split('\t')
		name = line_l[2].split()[0]
		all_cosmid[ name.split('.')[ 0 ]   ][name] = ''
		if name in no_unique:
			continue
		e_val = float( line_l[12] )
		if name not in all_blast_detail and e_val <= 1e-100:
			frame = line_l[18]
			if frame =='-1':
				tag = ' [RC] '
			else:
				tag = ' [] '
				
			align_length = abs( int( line_l[16] ) - int(  line_l[15] )  )
			all_blast_detail[ name ] =[line_l[6].split()[0], '#%s(0)'%( name )+tag +'%s bases, 00000000 checksum. '%( align_length  ) +'{%s %s} <1 %s>'%( line_l[15], line_l[16] ,  align_length   ),
			                            
			             sorted(   [int( line_l[16] ) , int(  line_l[15] )  ] ),
			             e_val
			             ]
		elif e_val >= 1e-100:
			continue
		else:
			if e_val==0:
				no_unique[ name ] = ''
			elif	all_blast_detail[ name  ][-1]==0:
				continue
					
			elif e_val/all_blast_detail[ name  ][-1] <=100:
				no_unique[ name ] = ''
				
				
	seq_blast = Ddict()

	for each_name in all_blast_detail:
		if each_name in no_unique:
			continue
		else:
			seq_blast[  all_blast_detail[ each_name  ][0]  ][  each_name ] = [   all_blast_detail[ each_name  ][1] ,all_blast_detail[ each_name  ][-2]
			]
			
	all_reference = {}


for each_seq in all_seq:
	#print( each_seq)
	if each_seq in seq_blast:
		OUTPUT.write( '##'+each_seq+' %s %sbases'%( len( seq_blast[ each_seq  ]  )  , 
		                                             len(  all_seq[ each_seq  ] )
		                                             )  
		             +', 00000000 checksum.\n' 
		              )
		OUTPUT.write( all_seq[ each_seq  ]+'\n'  )
		sequence = all_seq[ each_seq ]
		for each_read in seq_blast[ each_seq  ]:
			result = seq_blast[ each_seq  ][ each_read ]
			OUTPUT.write( result[0]+'\n'     )
			align = sequence[ result[-1][0] :result[-1][-1]      ]+'\n'
			OUTPUT.write( align )
	else:
		OTHER.write( '>'+each_seq+'\n' + all_seq[  each_seq]+'\n' )

PAIR.write(  'library	fosmid	2000	8000000\n'  )
REFBLAST.next()
for line in REFBLAST:
	line_l = line.split('\t')
	name = line_l[2].split()[0]
	if name in all_reference:
		no_unique[ name ] = ''
		continue
	all_reference[ name ] = ''
for each_name in all_cosmid:
	if len( all_cosmid[ each_name ]  )!=2:

		continue
	status = 'yes'
	for each_cos in all_cosmid[ each_name ]:
		if each_cos in no_unique:
			print( each_cos )
			status='no'
	if status =='yes':
		PAIR.write( '\t'.join( all_cosmid[ each_name ]   ) +'\tfosmid\n'     )
