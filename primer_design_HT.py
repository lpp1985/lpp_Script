#!/usr/bin/python
#Author=LPP
from lpp import *
import subprocess, string, shutil,multiprocessing
def check_path( path ):
	if os.path.exists(path):
		shutil.rmtree(path)
	os.mkdir( path )
	
def cout_Tm_f( gap_start ):
	global all_coords, s1
	if len([x for x in xrange( gap_start - 20  ,gap_start+1 ) if x in all_coords ]   )<20:
		return False
	else:
		seq = s1[ gap_start - 20  :   gap_start  ]
		if 'G' in seq:
			G_perc = seq.count( 'G' )
		else:
			G_perc = 0.0
		if 'C' in seq:
			C_perc = seq.count( 'C' )
		else:
			C_perc = 0.0
		G_C_P = (G_perc+C_perc)/float( len(seq) )
		if G_C_P <0.4:
			return False
		Tm = (G_perc+C_perc)*4+ 2*(  20 - G_perc+C_perc )
		return seq, str(Tm), str( G_C_P )
def cout_Tm_r( gap_stop ):
	global all_coords, s1

	if len( [x for x in xrange( gap_stop   ,gap_stop+20 ) if x in all_coords ]   )<20:
		return False
	else:
		seq = complement( s1[  gap_stop : gap_stop+20  ] )
		if 'G' in seq:
			G_perc = seq.count( 'G' )
		else:
			G_perc = 0.0
		if 'C' in seq:
			C_perc = seq.count( 'C' )
		else:
			C_perc = 0.0
		G_C_P = (G_perc+C_perc)/float( len(seq) )
		if G_C_P <0.4:
			return False
		Tm=(G_perc+C_perc)*4+ 2*(  20 - G_perc+C_perc )
		return seq,str(Tm), str( G_C_P )
path = './tag/'
check_path( path )
primer_Len = 20
frame_L = 8
def complement(char):
	libary=string.maketrans('atcgATCG','tagcTAGC')
	return char[::-1].translate( libary )
RAW = fasta_check(  open(sys.argv[1],'rU')  )




print('''Bowtie Build is running .....'''),
subprocess.call( ['bowtie-build' ,sys.argv[1], 'TOTAL'] ,stdout = subprocess.PIPE)
print('''ok''')



for t,s in RAW:
	title = re.search( '>(\S+)',t ).group(1)
	PRIMER = open( '%s.primer'%(  title ) ,'w')
	PRIMER.write('gap_ID	primer	site	length	Tm	GC%	length	gap_region\n')
	s1 = re.sub('\s+','',s)
	s2 = complement(  s1  )
	tag_hash = {}
	all_coords = {}
	for i in xrange(0,len( s1 )  ):
		all_coords[ i ] =''
	for i in xrange(frame_L,len(s1)+1):
		tag1 = s1[ i-frame_L:i  ]
		if not re.search('(N)',tag1,re.I):
			tag_hash[  tag1 ] = ''
		tag2 = s2[ i-frame_L:i  ]
		if not re.search('(N)',tag2,re.I):
			tag_hash[ tag2 ] = ''
	END= open( path+'%s.tag'%( title )  ,'w')
	i=0
	for tag in tag_hash:
		i=i+1
		END.write(  '>%s\n%s\n'%( i,tag  )  )
	END.close()
	pop = subprocess.Popen(  ['bowtie','-a','-n','0','TOTAL','-f',END.name ],stdout = subprocess.PIPE )
	end_hash = {}
	for line in pop.stdout:
		line_l = line.split('\t')
		number = line_l[0]
		if int( line_l[-2] ) !=0:
			if line_l[1] == '+':
				for k in xrange( int( line_l[-2] ) , int( line_l[-2] )+frame_L):
					if k in all_coords:
						del all_coords[ k ]
	i=0
	gap_tag = {}
	for each_Gap in re.finditer('(N+)',s1,re.I):
		i+=1
		[ gap_start, gap_stop ] = [ int(x) for x in each_Gap.span() ]
		gap_tag[ gap_start ] = [ title+'_gap%s'%( i ),  gap_stop  ]
		for j in  xrange(gap_start,gap_stop ):
			del all_coords[ j ]
	already={}
	for each_Gap in re.finditer('(N+)',s1,re.I):
		[ gap_start, gap_stop ] = [ int(x) for x in each_Gap.span() ]
		gap_start = gap_start - 400
		gap_stop  = gap_stop + 400
		if gap_start in already:
			continue
		while not cout_Tm_f(  gap_start  ):
			gap_start-=1
		primer_f,tm_f , G_C_P_f = cout_Tm_f(  gap_start  )
		 
		product_start = gap_start - primer_Len
		while not cout_Tm_r(  gap_stop  ):
			gap_stop+=1
		primer_r,tm_r , G_C_P_r = cout_Tm_r(  gap_stop  )
		product_stop = gap_stop + primer_Len 
		product_size = '%s'%(product_stop - product_start)
		tag_contain=[]
		gap_region = []
		for i in xrange( product_start,  product_stop  ):
			if i in gap_tag:
				already[ i ] = ''
				tag_contain.append( gap_tag[i][0] )
				gap_region.append('%s-%s'%( i - product_start,gap_tag[i][1] - product_start ) )
				
		PRIMER.write(  '; '.join( tag_contain )+'\t'+primer_f+'\t'+str(product_start)+'\t'+str( primer_Len )+'\t'+ tm_f+'\t'+ G_C_P_f+'\t'+product_size+'\t'+'; '.join( gap_region )+'\n\t')
		PRIMER.write(primer_r+'\t'+str(gap_stop)+'\t'+str( primer_Len )+'\t'+tm_r+'\t'+G_C_P_r+'\t'+product_size+'\t'+'; '.join( gap_region )+'\n')
