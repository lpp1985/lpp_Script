#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/4/24

from lpp import *
from optparse import OptionParser

status = Ddict()
if __name__=='__main__':
	usage = '''usage: python2.7 %prog [options] '''
	parser = OptionParser(usage =usage )



	parser.add_option("-e", "--EXONERATE", action="store",
		              dest="exo",
		              help="Exonerate_Result")

	parser.add_option("-o", "--out", action="store",
		              dest="out",
		              help="The output path  you want!!")
	parser.add_option("-r", "--raw", action="store",
		              dest="nucl",
		              help="query's nucelotide sequence!!")	
	(options, args) = parser.parse_args()
	EXO =block_reading( open( options.exo,'rU'  ),tag='C4 Alignment:'   ) 
	#MUM = open( options.mummer,'rU'  )
	EXO.next()
	nucleo= {}
	NUCLEO = fasta_check(  open( options.nucl,'rU'   )    )
	for t,s in NUCLEO:
		nucleo[ t[1:-1].split()[0].split(':')[-1]] = re.sub( '\s+','',s    )
	already = {}	
	all_has = Ddict()
	def exonerate_parse( each_b  ):
		#global status
		try:
			all_align = re.split( '\s+Target range\: [^\n]+\n\n',each_b  )[-1]
			all_align = re.split( '\s+cigar.+',all_align,re.DOTALL )[0]
		except:
			print( each_b )
		align_block = all_align.split('''\n\n''')
	#	try:
	#		direct = re.search( '([+-])	\.	\n',each_b ).group(1)
	#	except:
	#		continue
		frame_length = re.findall( '(\#+)',all_align  )
		start = 0
		seq_end = int(re.search(  'Query range\: \d+ \-\> (\d+)',each_b).group(1))
		for each_block in align_block:
			align_detail_block = each_block.split('\n')
			if start ==0:
				loca = re.search( '\w',align_detail_block[2]  ).span()[0]
				aa_raw = ''
				seq = ''
				align_tag = ''
				aa_new = ''
				query_start = int(re.search( '(\d+)\s+\:',align_detail_block[0] ).group(1) )
				query_start_aa = query_start
				if query_start!=1:
					query_start = query_start*3
				
			
			
			seq += re.sub( '[^ATCG\-]+','',  align_detail_block[-1]      )
			aa_raw += re.sub( '(?:^\s+\d+\s+\:\s+|\s+\:\s+\d+$)','', align_detail_block[0]   )
			aa_new += re.sub( '(?:^\s+|\s+$)','', align_detail_block[2]   )
			
			
			align_tag+= re.sub('(?:^\s{%s}|\n$)'%( loca ),'', align_detail_block[1])
			start = 1	

		if frame_length :
			all_frame = re.finditer( '(\#+)',align_tag  )
			for each_frame in all_frame:

				offset = each_frame.span()[0]
				
				status[ each_query  ][ 'Frame-shift'  ][ '.-->%s,%s'%( seq[ each_frame.span()[0]: each_frame.span()[1]  ]   ,query_start +offset )   ] = ''
				
		if 'Target Intron' in aa_raw:
			status[ each_query  ][ 'Target Intron'  ][ '-'   ] = ''
			#continue
		
		all_diff = re.finditer( '([^\|\#]+)',align_tag  )
		
		for each_diff in  all_diff:
			offset = each_diff.span()[0]
			
			diff_loca = query_start+each_diff.span()[0]
			
			if seq[:3] not in [ 'ATG','ATT','ATA','TTG','GTG'  ] or seq_end< len( seq[ each_diff.span()[0]: each_diff.span()[1]   ]  )/3:
				tag = 'Start-codon-break'
				if nucleo[ each_query.split()[0]  ][   each_diff.span()[0]: each_diff.span()[1]  ].upper() ==seq[ each_diff.span()[0]: each_diff.span()[1]   ].upper():
					continue
				status[ each_query  ][ tag  ][ '%s-->%s,%s'%( nucleo[ each_query.split()[0]  ][   each_diff.span()[0]: each_diff.span()[1]  ],seq[ each_diff.span()[0]: each_diff.span()[1]   ]   ,query_start +each_diff.span()[0]  )   ] ='-->'.join(  [aa_raw[  each_diff.span()[0]: each_diff.span()[1]     ] , aa_new[ each_diff.span()[0]: each_diff.span()[1]   ]] )
				break
			else:
				tag = 'Non-Anonoymious'
				
			status[ each_query  ][ tag  ][ '%s-->%s,%s'%( nucleo[ each_query.split()[0]  ][   each_diff.span()[0]: each_diff.span()[1]  ],seq[ each_diff.span()[0]: each_diff.span()[1]   ]   ,query_start +each_diff.span()[0]  )   ] ='-->'.join(  [aa_raw[  each_diff.span()[0]: each_diff.span()[1]     ] , aa_new[ each_diff.span()[0]: each_diff.span()[1]   ]] ) +',%s'%( query_start_aa+ each_diff.span()[0]/3  )		
	for each_b in EXO:
		
		
		query_name = re.search( 'Query\:\s+(.+?)\n' ,each_b   ).group(1)
		score = int( re.search( 'Raw score\:\s+(\d+)',each_b  ).group(1) )
		all_has[ query_name  ][score] = each_b
	for each_query in all_has:
		each_b = all_has[each_query][ sorted(  all_has[each_query]   )[-1]        ]
		exonerate_parse(each_b)
	
				
	END_SILENCE = open( options.out+'.silence' ,'w' )
	END_SILENCE.write( '\t'.join(['Gene_ID','Gene_annotation','Category', 'Why' ]) +'\n'  )
	END_MUTATION = open( options.out+'.mutation' ,'w' )
	END_MUTATION.write( '\t'.join(['Gene_ID','Gene_annotation', 'Non-Anonoymious',"Non-Anonoymious situation"] )+'\n'  )

	for each_gene in sorted(status):
		for key in status[each_gene]:
			if key =='Non-Anonoymious':
				for key1 in status[each_gene][ key ]: 
					END_MUTATION.write( '\t'.join([ each_gene.split()[0] ,  each_gene , key1,status[each_gene][key][key1]  ]    )+'\n' )
			else:
				for key1 in status[each_gene][ key ]: 
					END_SILENCE.write( '\t'.join([ each_gene.split()[0] ,  each_gene ,key ,key1  ]  )+'\n' )
			
