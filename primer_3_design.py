#!/usr/bin/python
#Author=LPP
from lpp import *
import shutil
import subprocess
import multiprocessing
def Tm_calc( seq ):
	gc_ctotain = len( re.findall( '([G|C])',seq ) ) 
	tm = 64.9+41*( gc_ctotain -16.4)/len(seq)
	if len(seq)>20:
		tm = tm+3
	return '%s'%( tm )
def check_path( path ):
	if os.path.exists(path):
		shutil.rmtree(path)
	os.makedirs( path )
def get_primer( (t,s) ):
	title =re.search( '>(\S+)',t ).group(1)
	s = re.sub('\s+','',s)
	all_gap_loc = re.finditer( '(N+)',s )
	all_n = re.findall( '(N+)',s )
	print( all_gap_loc )
	if len( all_n )<1:
		print(  '''$  %s is no need to design!!'''%(  title  )  )
		return ''
	print(  '''$  %s is  designing!!'''%(  title  )  )
	path = root_path+'/primer_input1/%s/'%( title )
	check_path(path)
	primer_end = path+title+'.primer'
	primer_raw = path+title+'.primer3'
	primer_log = path+title+'.log'

	s= re.sub('\s+','',s)
	all_gap = []


	END = open( primer_raw,'w' )
	LOG = open( primer_log,'w' )
	LOG.write('ID\tFrom\tTo\tProduct_Size\n')
	for each_gap in all_gap_loc:

		gap_start,gap_stop = each_gap.span()
		length = gap_stop - gap_start
	
		END.write( '''SEQUENCE_ID=%s
SEQUENCE_TEMPLATE=%s
SEQUENCE_TARGET=%s
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=17
PRIMER_MAX_SIZE=23
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=400-3000
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_MAX_TM=70
PRIMER_MIN_TM=30
PRIMER_OPT_TM=60
PRIMER_MAX_DIFF_TM=2
PRIMER_MAX_GC=65
PRIMER_MIN_GC=52
PRIMER_MAX_END_STABILITY=8
PRIMER_SELF_ANY=8
PRIMER_SELF_END=3
PRIMER_GC_CLAMP=0
PRIMER_PAIR_MAX_COMPL_ANY=4
PRIMER_PAIR_MAX_COMPL_END=3
=
'''%( title,s,'%s,%s'%( gap_start,length  )  ) )
		LOG.write( '''%s\t%s\t%s\t%s\n'''%( title,gap_start,gap_stop,length ) )
	os.chdir( path )
	subprocess.call(  ['primer3_core',primer_raw,'-output',title+'.primer'],stdout= subprocess.PIPE)
	RAW = block_reading( open(primer_end,'rU'),tag='^=' )
	no_hash = {}
	CANDIT = open( primer_log.replace('.log','.cand'),'w' )
	NO = open(  primer_log.replace('.log','.no'),'w' )
	CANDIT.write('gap	id\tprimer	site	length	Tm	GC%	length\tgap_from\tgap_to\tgap_region\n')
	CANDIT_FASTA = open(primer_log.replace('.log','.candfasta'),'w')
	i=0
	for e_b in RAW:
		i+=1
		gap_id = re.search( 'SEQUENCE_ID=(\S+)',e_b ).group(1)+'_gap%s'%( i )
		if re.search( '(PRIMER_LEFT_\d+_SEQUENCE)',e_b ):
			CANDIT.write( '%s'%(  gap_id  ) )
			gap_from , gap_length = re.search(  'SEQUENCE_TARGET=(\d+),(\d+)' ,e_b ).group( 1 , 2 )
			all_primer_block = re.split('PRIMER_PAIR_NUM_RETURNED=\d+',e_b)[-1]
			primers_data_list = re.split('PRIMER_PAIR_\d_PENALTY=\S+\n',all_primer_block)[1:]
			for each_primer_data in primers_data_list:
				id__left = re.search('PRIMER_LEFT_(\d+)_SEQUENCE=(\w+)',each_primer_data)
				p_id,leftseq = id__left.group(1),id__left.group(2)
				rightseq = re.search('PRIMER_RIGHT_\d+_SEQUENCE=(\w+)',each_primer_data).group(1)
				CANDIT_FASTA.write( '>%s_primer%s__left\n%s\n'%( gap_id ,p_id,leftseq ) )
				CANDIT_FASTA.write( '>%s_primer%s__right\n%s\n'%( gap_id ,p_id,rightseq ) )
				left_data= re.search(  'PRIMER_LEFT_\d+=(\d+),(\d+)'  , each_primer_data )
				left_site,left_length = left_data.group(1),left_data.group(2)
				right_data= re.search(  'PRIMER_RIGHT_\d+=(\d+),(\d+)'  , each_primer_data )
				left_site,left_length = left_data.group(1),left_data.group(2)
				right_site,right_length = right_data.group(1),right_data.group(2)
				tm_gc = re.search('''PRIMER_LEFT_\d+_TM=(\S+)
PRIMER_RIGHT_\d+_TM=(\S+)
PRIMER_LEFT_\d+_GC_PERCENT=(\S+)
PRIMER_RIGHT_\d+_GC_PERCENT=(\S+)''',each_primer_data )
				tm_left ,tm_right,gc_left,gc_right = tm_gc.group(1,2,3,4)
				tm_left = Tm_calc(  leftseq  )
				tm_right = Tm_calc(  rightseq  )
				product_size = re.search('PRIMER_PAIR_\d+_PRODUCT_SIZE=(\d+)',each_primer_data).group(1)
				
				gap_data = '%s\t%s\t%s--%s\n'%( gap_from,int(gap_from)+int( gap_length ), int(gap_from)- int(left_site),int(gap_from)- int(left_site)+int( gap_length ) )
				
				output_gap1 = '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(   p_id,leftseq,left_site,left_length,tm_left,gc_left,product_size, gap_data)
				output_gap2 = '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(   p_id,rightseq,right_site,right_length,tm_right,gc_right,product_size, gap_data )
				CANDIT.write( output_gap1 + output_gap2 )
		else:
			no_hash[ gap_id ] = ''
	LOG = open(glob.glob('*.log')[0])
	for line in LOG:
		g_id = line.split('\t')[0]
		if g_id in no_hash:
			NO.write(line)
	pop = subprocess.Popen(  ['bowtie',root_path+'TOTAL','-f',primer_log.replace('.log','.candfasta') ],stdout= subprocess.PIPE)
	already= {}
	no_dup = {}
	ERROR = open( path + title+'.error','w')
	for line in pop.stdout:
		r_id = line.split()[0]
		if r_id not in already:
			no_dup [ r_id ] = ''
		if r_id in no_dup [ r_id ]:
			del no_dup [ r_id ]
			ERROR.write( line )
		already[ r_id ] = ''
	
	print('!  %s designing process is finished'%( title )  )


root_path = os.path.abspath('./')+'/'
ALL_PRIMER = open('all_primer.list','w')
BIG_GAP = open('BIG_primer.list','w')
ALL_CAND = open('all_primer.cand','w')

BIG_CAND = open('BIG_primer.cand','w')


print('''Bowtie Build is running....''')
subprocess.call( ['bowtie-build' ,sys.argv[1], 'TOTAL'],stdout = subprocess.PIPE )

pool = multiprocessing.Pool( processes=20 ) 
RAW = fasta_check(  open(sys.argv[2],'rU')  )
print( 'Start to design Primer......' )

pool.map( get_primer,  RAW )
for roots,dirs,files in os.walk(root_path+'/primer_input1/'):
	for e_f in files:
		if e_f.endswith('.cand'):
			RAW = open( roots+'/'+e_f )
			RAW.next()
			for line in RAW:
				line_l = line.split('\t')
				if line_l[0]:
					seq_id = line_l[0]
					tag = 'F'
				else:
					tag = 'R'
				if line_l[1]=='0':
					ALL_PRIMER.write(  seq_id+' '+tag +'\t'+ line_l[2]+'\n')
					start,stop = re.search( '(\d+)\-\-(\d+)' ,line).group(1), re.search( '(\d+)\-\-(\d+)' ,line).group(2)
					gap_len = int(stop)-int(start)
					if gap_len>1000:
						BIG_GAP.write( seq_id+' '+tag +'\t'+ line_l[2]+'\n' )
						BIG_CAND.write(line)
					ALL_CAND.write(line)
		
