#!/usr/bin/python
#Author=LPP
from lpp import *
import shutil
import subprocess
import multiprocessing
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly predict primer among different template!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-I", "--MIN", action="store", 
                  dest="min",
                  type='string',
                  default = '400',
                  help="the minium size of product [ 400 ]")

parser.add_option("-X", "--MAX", action="store", 
                  dest="max",
                  type='string',  
                  default = '3000',
                  help="the maximum size of product [ 3000  ]")
parser.add_option("-R", "--reference", action="store", 
                  dest="reference",
                  type='string',  
                  help="the reference sequence for primer specificity testing by bowtie")

parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  type='string',  
                  help="the sequence you want to design primer")

parser.add_option("-U", "--UPPER", action="store", 
                  dest="upper",
                  type='string', 
                  default = '60',
                  help="the maximum Tm temperature allowed")

parser.add_option("-D", "--DOWN", action="store", 
                  dest="down",
                  type='string',
                  default = '40',
                  help="the miimum Tm temperature allowed")

parser.add_option("-O", "--OPTIMUM", action="store", 
                  dest="optimum",
                  type='string',
                  default = '50',
                  help="the optimum Tm temperature allowed")

parser.add_option("-d", "--DIFF", action="store", 
                  dest="diff",
                  type='string',
                  default = '2',
                  help="the difference of Tm between a pair of primer" )

parser.add_option("-G", "--G_UP", action="store", 
                  dest="g_upper",
                  type='string',
                  default = '65',
                  help="the maximum of GC% of primer" )

parser.add_option("-g", "--G_down", action="store", 
                  dest="g_down",
                  type='string',
                  default = '30',
                  help="the minimum of GC% of primer" )

(options, args) = parser.parse_args() 

'''Start to design'''


product_min = options.min
product_max = options.max

reference = options.reference
template = options.input

Tm_upper = options.upper
Tm_dwon = options.down
Tm_optimum = options.optimum
Tm_diff = options.diff

gc_upper = options.g_upper
gc_down = options.g_down


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
	all_gap_loc = re.finditer( '([n|N]+)',s )
	all_n = re.findall( '([n|N]+)',s )
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
	i=0
	for each_gap in all_gap_loc:
		i+=1
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
PRIMER_PRODUCT_SIZE_RANGE=%s-%s
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_MAX_TM=%s
PRIMER_MIN_TM=%s
PRIMER_OPT_TM=%s
PRIMER_MAX_DIFF_TM=%s
PRIMER_MAX_GC=%s
PRIMER_MIN_GC=%s
PRIMER_MAX_END_STABILITY=8
PRIMER_SELF_ANY=8
PRIMER_SELF_END=3
PRIMER_GC_CLAMP=0
PRIMER_PAIR_MAX_COMPL_ANY=4
PRIMER_PAIR_MAX_COMPL_END=3
=
'''%( title,s,'%s,%s'%( gap_start,length  ),product_min,product_max,Tm_upper,Tm_dwon,Tm_optimum, Tm_diff, gc_upper,gc_down  ) )
		LOG.write( '''%s\t%s\t%s\t%s\n'''%( title+'_%s'%(i),gap_start,gap_stop,length ) )
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
		gap_id = re.search( 'SEQUENCE_ID=(\S+)',e_b ).group(1)+'_%s'%( i )
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
	print( no_hash )
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
ALL_CAND = open('all_primer.cand','w')



print('''Bowtie Build is running....''')
subprocess.call( ['bowtie-build' ,reference, 'TOTAL'],stdout = subprocess.PIPE )

pool = multiprocessing.Pool( processes=20 ) 
RAW = fasta_check(  open(  template,'rU'  )  )
print( 'Start to design Primer......' )
NO = open( 'Primer3.fail','w' )
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
				if line_l[1]=='0' :
					ALL_PRIMER.write(  seq_id+tag +'\t'+ line_l[2]+'\n')
					start,stop = re.search( '(\d+)\-\-(\d+)' ,line).group(1), re.search( '(\d+)\-\-(\d+)' ,line).group(2)
					gap_len = int(stop)-int(start)
					line_l[1] = seq_id +tag
					line = '\t'.join( line_l )
					ALL_CAND.write(line)
		elif e_f.endswith( '.no' ):
			RAW_NO = open( roots+'/'+e_f )
			for line in RAW_NO:
				NO.write( line.split()[0]+'\n' )
		