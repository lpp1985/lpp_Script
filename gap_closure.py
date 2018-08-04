#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/13
from lpp import *
from optparse import OptionParser 
usage = '''usage: python2.7 %prog -s sequence_result( 454,or 3730  ) -l gap_location -b blast_result'''
parser = OptionParser(usage =usage )
ok_closure = {}
parser.add_option("-s", "--SEQUENCE", action="store", 
                  dest="seq", 
                  help="Sequencing result for gap closure")
parser.add_option("-l", "--LOCATION", action="store", 
                  dest="location", 
                  help="gap location file")
parser.add_option("-b", "--BLAST", action="store", 
                  dest="blast", 
                  help="Blast result(  sequencing result( query ) VS scaffold result( subject ) )")
(options, args) = parser.parse_args() 
seq_file   = options.seq

loca_file  = options.location

blast_file  = options.blast

SEQ = fasta_check( open( seq_file,'rU'  ) )
all_seq = {}
for t,s in SEQ:
	t = t[1:-1]
	s = re.sub('\s+','',s)
	all_seq[ t ] = s
	


RAW = open( loca_file,'rU'  )
all_location = Ddict()
for line in RAW:
	line_l = line[1:-1].upper().split('\t')
	
	all_location[ line_l[0]  ][ int(  line_l[1] )   ] = ''
	all_location[ line_l[0]  ][ int(  line_l[-1] )   ] = ''
BLAST = open( blast_file ,'rU' )
line_number = int( re.search('^(\d+)',  os.popen( 'wc -l %s'%(  blast_file   ) ).read() ).group(1) )
blast_location = Ddict()
blast_frame = Ddict()
already_close = {}
i=0
LOG = open( 'closure_situation.log' ,'w' )
FASTA = open( 'alread_closure.fasta','w'  )
BLAST.next()
line = BLAST.next()
line_l = line[:-1].split('\t')
gap_id = re.search( '(s\d+\-\d+)',line ).group(1)
i=i+1
if 'scaffold'  in line:
	hit_scaffold = re.search( 'scaffold(\d+)',line ).group(1)
	gap_source = re.search( 'S(\d+)',line  ).group(1)
	
	already_close[ line_l[2]  ] = ''
	if gap_source == hit_scaffold:
		blast_location[ gap_id ][  int( line_l[ 15 ] )  ]= int( line_l[ 13 ] )  
		blast_location[ gap_id ][  int( line_l[ 16 ] )  ]= int( line_l[ 14 ] ) 
		blast_frame[  int( line_l[ 18 ] ) ][ int( line_l[ 15 ] )  ] = ''
		blast_frame[  int( line_l[ 18 ] )  ][  int( line_l[ 16 ] )  ] = ''
		query_length = int( line_l[3]  )
		query_id = line_l[2]

for line in BLAST:

	i=i+1
	
	if 'scaffold' not in line:
		tqg=1
		continue
	hit_scaffold = re.search( 'scaffold(\d+)',line ).group(1)
	gap_source = re.search( 'S(\d+)',line  ).group(1)
	line_l = line[:-1].split('\t')
	

	if i == line_number or line_l[2] not in already_close :
		for each_gap in blast_location:
			
			if each_gap in ok_closure:
				continue
			[gap_from,gap_end] = sorted(  all_location[ each_gap ].keys()   )
			gap_compare = {}
			for each_frame in blast_frame:
				distance_from = sorted( [  each_location - gap_from    for each_location in  blast_frame[ each_frame ]       ]  ,key = lambda x: abs(x) )[0]
				distance_end = sorted( [  each_location - gap_end    for each_location in  blast_frame[ each_frame ]       ]  ,key = lambda x: abs(x))[0]
					
				output_from = blast_location[ each_gap ][ gap_from + distance_from ]  - ( distance_from * each_frame ) 
				output_end  = blast_location[ each_gap ][ gap_end + distance_end ]- ( distance_end * each_frame )
				[print_from,print_end] = sorted([ output_from, output_end   ] )
				print_from+=1; print_end+=1
				gap_compare[  abs( distance_from ) + abs( distance_end ) ] = [ output_from, output_end , query_id]
				if output_from > output_end:
					output_frame = '-'
					[output_from,output_end] = sorted([ output_from, output_end   ] )
					if output_from <0 and output_end > query_length:
						status = 'no gap found'
					elif output_from <0:
						status = 'only found gap tail'
					elif output_end > query_length:
						status = 'only found gap head'
					else:
						status = 'gap_found'
						ok_closure[  each_gap ] = ''
						sequence = complement( all_seq[ query_id ][ output_from: output_end  ] )
						FASTA.write( '>'+each_gap+'\n'+sequence+'\n'  )
					LOG.write( '%s\t%s\t%s\t%s\t%s\n'%( query_id,output_frame ,output_from,output_end ,status )    )
				else:
					
					output_frame = '+'
					if print_from <0 and print_end > query_length:
						status = 'no gap found'
					elif print_from <0:
						status = 'only found gap tail'
					elif print_end > query_length:
						status = 'only found gap head'
					else:
						status = 'gap_found'
						sequence =  all_seq[ query_id ][ print_from: print_end  ] 
						ok_closure[  each_gap ] = ''
						FASTA.write( '>'+each_gap+'\n' +sequence+'\n' )
					LOG.write( '%s\t%s\t%s\t%s\t%s\n'%( query_id,output_frame ,print_from,print_end ,status )    )
		blast_location = Ddict()
		blast_frame = Ddict()
		gap_id = re.search( '(S\d+\-\d+)',line ).group(1)
		already_close[ line_l[2]  ] = ''
		if gap_source == hit_scaffold:
			blast_location[ gap_id ][  int( line_l[ 15 ] )  ]= int( line_l[ 13 ] )  
			blast_location[ gap_id ][  int( line_l[ 16 ] )  ]= int( line_l[ 14 ] ) 
			blast_frame[  int( line_l[ 18 ] ) ][ int( line_l[ 15 ] )  ] = ''
			blast_frame[  int( line_l[ 18 ] )  ][  int( line_l[ 16 ] )  ] = ''
			query_length = int( line_l[3]  )
			query_id = line_l[2]
	else:

		gap_id = re.search( '(S\d+\-\d+)',line ).group(1)
		
		
		if gap_source == hit_scaffold:
			blast_location[ gap_id ][  int( line_l[ 15 ] )  ]= int( line_l[ 13 ] )  
			blast_location[ gap_id ][  int( line_l[ 16 ] )  ]= int( line_l[ 14 ] ) 
			blast_frame[  int( line_l[ 18 ] ) ][ int( line_l[ 15 ] )  ] = ''
			blast_frame[  int( line_l[ 18 ] )  ][  int( line_l[ 16 ] )  ] = ''
			query_length = int( line_l[3]  )
			query_id = line_l[2]
else:
	for each_gap in blast_location:
			
		if each_gap in ok_closure:
			continue
		[gap_from,gap_end] = sorted(  all_location[ each_gap ].keys()   )
		gap_compare = {}
		for each_frame in blast_frame:
			distance_from = sorted( [  each_location - gap_from    for each_location in  blast_frame[ each_frame ]       ]  ,key = lambda x: abs(x) )[0]
			distance_end = sorted( [  each_location - gap_end    for each_location in  blast_frame[ each_frame ]       ]  ,key = lambda x: abs(x))[0]
				
			output_from = blast_location[ each_gap ][ gap_from + distance_from ]  - ( distance_from * each_frame ) -1
			output_end  = blast_location[ each_gap ][ gap_end + distance_end ]- ( distance_end * each_frame ) -1
			gap_compare[  abs( distance_from ) + abs( distance_end ) ] = [ output_from, output_end , query_id]
			if output_from > output_end:
				output_frame = '-'
				[output_from , output_end]  = sorted(  [output_from , output_end]   )
				if output_from <0 and output_end > query_length:
					status = 'no gap found'
				elif output_from <0:
					status = 'only found gap tail'
				elif output_end > query_length:
					status = 'only found gap head'
				else:
					status = 'gap_found'
					ok_closure[  each_gap ] = ''
					sequence = complement( all_seq[ query_id ][ output_from: output_end  ] )
					FASTA.write( '>'+each_gap+'\n'+sequence+'\n'  )
				LOG.write( '%s\t%s\t%s\t%s\t%s\n'%( query_id,output_frame ,output_from,output_end ,status )    )
			else:
				
				output_frame = '+'
				if output_from <0 and output_end > query_length:
					status = 'no gap found'
				elif output_from <0:
					status = 'only found gap tail'
				elif output_end > query_length:
					status = 'only found gap head'
				else:
					status = 'gap_found'
					sequence =  all_seq[ query_id ][ output_from+1: output_end+1  ] 
					ok_closure[  each_gap ] = ''
					FASTA.write( '>'+each_gap+'\n' +sequence+'\n' )
				LOG.write( '%s\t%s\t%s\t%s\t%s\n'%( query_id,output_frame ,output_from,output_end ,status )    )
				
FASTA = open(  FASTA.name,'rU'  )

already = dict( re.findall( '>(\S+)()',FASTA.read()   ) )
TOTAL = open(  RAW.name,'rU'  )
total = dict( re.findall( '>(\S+)()',TOTAL.read().upper()   ) )
error = {}
for key1 in total :
	if key1 not in already:
		error[key1] = ''
print(  len( error )  )
ERROR = open('error.log','w')
statc = Ddict()
all_stac = {}
blank = {}
for each_key in error:
	if not os.path.isdir(  '../trim/%s/'%( each_key )  ):
		blank[ each_key  ] = each_key+'\t'+'-\t-\t-\t-\tNO sequencing result!!!\n'
	for each_file in glob.glob( '../trim/%s/*.fasta'%( each_key ) ):
		if '(' in each_file:
			name = re.search( '([^-]+)\)',each_file ).group(1)
			name = re.sub( '^\d+','',name )
			statc[ each_key ][ name ] = ''
			all_stac[ name ] = ''
ERROR.write(  'ID'  )
for key1 in sorted( all_stac ):
	ERROR.write( '\t'+ key1)
ERROR.write( '\n' )
for key1 in statc:
	ERROR.write(  key1 )
	error_log = 'not specific'
	tag=[]
	if 'F' in statc[ key1 ] or 'R' in statc[ key1 ]:
		if 'F' not in statc[ key1 ] or 'R' not in statc[ key1 ]:
			tag += [x for x in [ 'F','R' ] if x not in statc[ key1 ] ]

	if 'F2' in statc[ key1 ] or 'R2' in statc[ key1 ]:
		if 'F2' not in statc[ key1 ] or 'R2' not in statc[ key1 ]:
			tag += [x for x in [ 'F2','R2','F','R' ] if x not in statc[ key1 ] ]
	if not tag:
		data = re.sub('\s+','',open( '../trim/%s/%s.fasta.contigs'%( key1, key1 ),'rU'   ).read())
		if '>' in data:
			error_log = 'not specific!!'
		else:
			error_log = 'need enlongation!!'
	else:
		error_log = '%s is error!!'%(' '.join(tag) )
	for each_key in sorted( all_stac ):
		if each_key in statc[  key1 ]:
			ERROR.write(  '\t*' )
		else:
			ERROR.write(  '\t-' )
	ERROR.write('\t'+error_log)
	ERROR.write('\n')
for key1 in blank:
	ERROR.write(blank[  key1 ])