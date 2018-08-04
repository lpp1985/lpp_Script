#!/usr/bin/python
#Author=LPP
from lpp import *
import shutil
import subprocess
import multiprocessing
from optparse import OptionParser
from email.mime.text import MIMEText
import cgi,glob,subprocess,time,smtplib,shutil
def send_mail( mail_address, job_title , download_link  ):
	toaddrs  = mail_address
	fromaddr = 'gasscine@gmail.com'
	msg = MIMEText('''Dear user
	    Your job named %s has been finished. You could download the result from %s 
	                                                                              Best wishes and enjoy your work!!!'''%( job_title, download_link    )  )
	msg[ 'Subject' ] = 'Your Misson in CHGB Blast has completed'
	msg[ 'From' ] = 'from@gmail.com'
	msg[ 'To' ] = mail_address

	# Credentials (if needed)
	username = 'gasscine'
	password = '28324086'

	# The actual mail send
	server = smtplib.SMTP('smtp.gmail.com')
	server.starttls()
	server.login(username,password)
	server.sendmail(fromaddr, toaddrs, msg.as_string())
	server.quit()


usage = '''usage: python %prog [options] 

It can automaticly predict primer among different template!!'''
parser = OptionParser(usage =usage ) 
parser.add_option("-I", "--MIN", action="store", 
                  dest="min",
                  type='string',
                  default = '400',
                  help="the minium size of product [ 400 ]")

parser.add_option("-o", "--output", action="store", 
                  dest="output",
                  type='string',
                  help="the output_path")

parser.add_option("-m", "--MAIL", action="store", 
                  dest="mail",
                  type='string',
                  help="the mail address ")
parser.add_option("-t", "--JOB_TITLE", action="store", 
                  dest="job_title",
                  type='string',
                  help="the job_title ")

parser.add_option("-X", "--MAX", action="store", 
                  dest="max",
                  type='string',  
                  default = '3000',
                  help="the maximum size of product [ 3000  ]")
parser.add_option("-r", "--REFERENCE", action="store", 
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
                  help="the maximum Tm temperature allowed[  60  ]")

parser.add_option("-D", "--DOWN", action="store", 
                  dest="down",
                  type='string',
                  default = '40',
                  help="the miimum Tm temperature allowed[  40  ]")

parser.add_option("-O", "--OPTIMUM", action="store", 
                  dest="optimum",
                  type='string',
                  default = '50',
                  help="the optimum Tm temperature allowed[  50  ]")

parser.add_option("-d", "--DIFF", action="store", 
                  dest="diff",
                  type='string',
                  default = '2',
                  help="the difference of Tm between a pair of primer [  2 ]" )

parser.add_option("-G", "--G_UP", action="store", 
                  dest="g_upper",
                  type='string',
                  default = '65',
                  help="the maximum of GC% of primer [  65  ]" )

parser.add_option("-g", "--G_DOWN", action="store", 
                  dest="g_down",
                  type='string',
                  default = '30',
                  help="the minimum of GC% of primer [  30  ]" )
parser.add_option("-v", "--OVERLAP", action="store", 
                  dest="overlap",
                  type='string',
                  default = '100',
                  help="the minimum size or overlap between product and Template [  100  ]" )
(options, args) = parser.parse_args() 

'''Start to design'''

overlap = int( options.overlap )
product_min = options.min
product_max = options.max

mail_address = options.mail
output_path = options.output
job_title = options.job_title

reference = options.reference
template = options.input

Tm_upper = options.upper
Tm_dwon = options.down
Tm_optimum = options.optimum
Tm_diff = options.diff

gc_upper = options.g_upper
gc_down = options.g_down

not_design = {}
def Tm_calc( seq ):
	gc_ctotain = len( re.findall( '([G|C])',seq ) ) 
	tm = 64.9+41*( gc_ctotain -16.4)/len(seq)
	if len(seq)>20:
		tm = tm+3
	return '%s'%( tm )
def check_path( path ):
	if not  os.path.exists(path):

		os.makedirs( path )
def get_primer( (t,s) ):
	title =re.search( '>(\S+)',t ).group(1)
	s = re.sub('\s+','',s)
	all_gap_loc = re.finditer( '([n|N]+)',s )
	all_n = re.findall( '([n|N]+)',s )
	if len( all_n )<1:
		#print(  '''$  %s is no need to design!!'''%(  title  )  )
		return ''
	#print(  '''$  %s is  designing!!'''%(  title  )  )
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
	os.system(  ' '.join(['primer3_core',primer_raw,'-output',title+'.primer' ,'>nul'] )  )
	#subprocess.call(  ['primer3_core',primer_raw,'-output',title+'.primer','>nul','2>&1'],stdout= subprocess.PIPE)
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

	for line in LOG:
		g_id = line.split('\t')[0]
		if g_id in no_hash:
			NO.write(line)
	ALL_CAND = open(  primer_log.replace('.log','.candfasta')  ,'rU' )
	already= {}
	no_dup = {}
	ERROR = open( path + title+'.error','w')
	if ALL_CAND.read().strip():
		pop = subprocess.Popen(  ['bowtie',root_path+'TOTAL','-p','50','-f','--quiet',primer_log.replace('.log','.candfasta') ],stdout= subprocess.PIPE)
		
		for line in pop.stdout:
	
			r_id = line.split()[0]
			if r_id not in already:
				no_dup [ r_id ] = ''
			if r_id in no_dup [ r_id ]:
				del no_dup [ r_id ]
				ERROR.write( line )
			already[ r_id ] = ''

	#print('!  %s designing process is finished'%( title )  )


root_path = os.path.abspath('./')+'/'
ALL_PRIMER = open('all_primer.list','w')
ALL_CAND = open('all_primer.cand','w')
ALL_CAND.write('gap_ID	primer	site	length	Tm	GC%	Product_length	gap_region\n')


#print('''Bowtie Build is running....''')
subprocess.call( ['bowtie-build' ,'-q',reference, 'TOTAL',],stdout = subprocess.PIPE )

pool = multiprocessing.Pool( processes=20 ) 
RAW = fasta_check(  open(  template,'rU'  )  )
#print( 'Start to design Primer......' )
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
				not_design[ line.split()[0] ] = ''



'''---Primer_desing.py------'''


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
RAW = fasta_check(  open(template,'rU')  )

FAIL = open( 'Fail_design_scan.log' ,'w'   )


for t,s in RAW:
	failed = {}
	title = re.search( '>(\S+)',t ).group(1)

	s1 = re.sub('\s+','',s)
	s2 = complement(  s1  )
	tag_hash = {}
	all_coords = {}
	for i in xrange(0,len( s1 )  ):
		all_coords[ i ] =''
	for i in xrange(frame_L,len(s1)+1):
		tag1 = s1[ i-frame_L:i  ]
		if not re.search('(n|N)',tag1,re.I):
			tag_hash[  tag1 ] = ''
		tag2 = s2[ i-frame_L:i  ]
		if not re.search('(n|N)',tag2,re.I):
			tag_hash[ tag2 ] = ''
	END= open( path+'%s.tag'%( title )  ,'w')
	i=0
	for tag in tag_hash:
		i=i+1
		END.write(  '>%s\n%s\n'%( i,tag  )  )
	END.close()
	pop = subprocess.Popen(  ['bowtie','-a','-p','50','--quiet','-n','0','TOTAL','-f',END.name ],stdout = subprocess.PIPE )
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
	all_gap = {}
	for each_Gap in re.finditer('([n|N]+)',s1):
		i+=1
		title_new = title+'_%s'%(i)
		if title_new not in not_design:
			continue

		[ gap_start_raw, gap_stop_raw ] = [ int(x) for x in each_Gap.span() ]

		gap_start = gap_start_raw -overlap
		gap_stop  = gap_stop_raw +overlap
		gap_tag[ gap_start_raw ] = [ title_new,  gap_stop_raw  ]
		all_gap[  title_new ] = ''
		for j in  xrange(gap_start_raw,gap_stop_raw ):
			if j in all_coords:
				del all_coords[ j ]
	already={}
	primer_total_log = {}
	k=0

	for each_Gap in re.finditer('([n|N]+)',s1,re.I):
		k+=1
		title_new = title+'_%s'%( k )

		if title_new not in not_design:

			continue
		[ gap_start_raw, gap_stop_raw ] = [ int(x) for x in each_Gap.span() ]
		gap_start = gap_start_raw -overlap
		gap_stop  = gap_stop_raw +overlap
		if gap_start in already:
			continue
		while not cout_Tm_f(  gap_start  ):
			gap_start-=1
			if gap_start ==1:
				break
		if gap_start ==1:
			continue
		primer_f,tm_f , G_C_P_f = cout_Tm_f(  gap_start  )

		product_start = gap_start - primer_Len
		while not cout_Tm_r(  gap_stop  ):
			gap_stop+=1
			if gap_stop == len(  s1  ):
				break
		if gap_stop == len(  s1  ):
			continue
		if gap_stop -gap_start>3000:
			continue
		primer_r,tm_r , G_C_P_r = cout_Tm_r(  gap_stop  )
		product_stop = gap_stop + primer_Len 
		product_size = '%s'%( product_stop - product_start )
		tag_contain=[]
		gap_region = []
		for i in xrange( product_start,  product_stop  ):
			if i in gap_tag:
				primer_total_log[ gap_tag[i][0] ] = ''
				already[ i ] = ''
				tag_contain.append( gap_tag[i][0] )
				gap_region.append('%s-%s'%( i - product_start,gap_tag[i][1] - product_start ) )

		ALL_CAND.write(  '; '.join( tag_contain )+'\t'+'; '.join(  [x+'F' for x in tag_contain ]  )  +'\t'+primer_f+'\t'+str(product_start)+'\t'+str( primer_Len )+'\t'+ tm_f+'\t'+ G_C_P_f+'\t'+product_size+'\t'+'; '.join( gap_region )+'\tHT\n\t')
		ALL_PRIMER.write(  '; '.join(  [x+'F' for x in tag_contain ]  )  +'\t'+primer_f+'\n'  )
		ALL_CAND.write('; '.join(  [x+'R' for x in tag_contain ]  )  +'\t'+primer_r+'\t'+str(gap_stop)+'\t'+str( primer_Len )+'\t'+tm_r+'\t'+G_C_P_r+'\t'+product_size+'\t'+'; '.join( gap_region )+'\tHT\n')
		ALL_PRIMER.write(  '; '.join(  [x+'R' for x in tag_contain ]  )  +'\t'+primer_f+'\n'  )
		already = {}
	for each_gap in not_design:
		if each_gap  in primer_total_log:
			already[ each_gap ] = ''

	for key1 in already:
		del not_design[ key1 ]



'''-----Final----'''
def cout_Tm_f_try( gap_start ):
	global all_coords, s1

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
def cout_Tm_r_try( gap_stop ):
	global all_coords, s1

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



for t,s in RAW:
	failed = {}
	title = re.search( '>(\S+)',t ).group(1)

	s1 = re.sub('\s+','',s)
	s2 = complement(  s1  )
	tag_hash = {}
	all_coords = {}

	i=0
	gap_tag = {}
	all_gap = {}
	for each_Gap in re.finditer('([n|N]+)',s1):
		i+=1
		title_new = title+'_%s'%(i)
		if title_new not in not_design:
			continue

		[ gap_start_raw, gap_stop_raw ] = [ int(x) for x in each_Gap.span() ]

		gap_start = gap_start_raw -overlap
		gap_stop  = gap_stop_raw +overlap
		gap_tag[ gap_start_raw ] = [ title_new,  gap_stop_raw  ]
		all_gap[  title_new ] = ''
		for j in  xrange(gap_start_raw,gap_stop_raw ):
			if j in all_coords:
				del all_coords[ j ]
	already={}
	primer_total_log = {}
	k=0

	for each_Gap in re.finditer('([n|N]+)',s1,re.I):
		k+=1
		title_new = title+'_%s'%( k )

		if title_new not in not_design:

			continue
		[ gap_start_raw, gap_stop_raw ] = [ int(x) for x in each_Gap.span() ]
		gap_start = gap_start_raw -overlap
		gap_stop  = gap_stop_raw +overlap
		if gap_start in already:
			continue

		primer_f,tm_f , G_C_P_f = cout_Tm_f_try(  gap_start  )

		product_start = gap_start - primer_Len

		primer_r,tm_r , G_C_P_r = cout_Tm_r_try(  gap_stop  )
		product_stop = gap_stop + primer_Len 
		product_size = '%s'%( product_stop - product_start )
		tag_contain=[]
		gap_region = []
		for i in xrange( product_start,  product_stop  ):
			if i in gap_tag:
				primer_total_log[ gap_tag[i][0] ] = ''
				already[ i ] = ''
				tag_contain.append( gap_tag[i][0] )
				gap_region.append('%s-%s'%( i - product_start,gap_tag[i][1] - product_start ) )

		ALL_CAND.write(  '; '.join( tag_contain )+'\t'+'; '.join(  [x+'F' for x in tag_contain ]  )  +'\t'+primer_f+'\t'+str(product_start)+'\t'+str( primer_Len )+'\t'+ tm_f+'\t'+ G_C_P_f+'\t'+product_size+'\t'+'; '.join( gap_region )+'\tTry\n\t')
		ALL_PRIMER.write(  '; '.join(  [x+'F' for x in tag_contain ]  )  +'\t'+primer_f+'\n'  )
		ALL_CAND.write('; '.join(  [x+'R' for x in tag_contain ]  )  +'\t'+primer_r+'\t'+str(gap_stop)+'\t'+str( primer_Len )+'\t'+tm_r+'\t'+G_C_P_r+'\t'+product_size+'\t'+'; '.join( gap_region )+'\tTry\n')
		ALL_PRIMER.write(  '; '.join(  [x+'R' for x in tag_contain ]  )  +'\t'+primer_f+'\n'  )

os.system( 'mv ./all_primer.* %s && rm %s -rf '%( output_path+'/','./' )   )
#os.system( 'mv ./all_primer.* %s'%( output_path+'/' )   )
send_mail( mail_address,job_title, output_path  )
