#!/usr/bin/env python
#coding:gb2312
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/6
from lpp import *
import sys, os, zipfile,shutil,glob,os
from optparse import OptionParser 
import subprocess
import time
usage = '''usage: python2.7 %prog -i input_path -o output_path'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  default = './', 
                  help="Input Path")
parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output", 
                  help="Output Path")
(options, args) = parser.parse_args() 
outpath = options.output

trace_outputpath = outpath+'trace/'
assembly_path = outpath+'assembly/'
cali_path = outpath+'caliber/'
trim_path = outpath+'trim/'
blast_path = outpath+'blast_check/'
intputpath = options.input
def combine_fasta(  ALL_SEQ ):
	all_seq = glob.glob( '*.fasta'  )
	for e_f in all_seq:
		data = open( e_f,'rU'  ).read()
		if not data.endswith(  '\n'  ):
			data = data+'\n'
		data = re.sub('\([^\(]+\)\.ab1\.seq','',data)
		ALL_SEQ.write( data  )  
def all_fasta(ALL_SEQ):
	RAW = open( cache_name,'rU')
	data = RAW.read()
	data = re.sub('\([^\(]+?\)\.ab1\.seq','',data)
	ALL_SEQ.write( data+'\n'  )
def check__trim(  path  ):
	def check_qual( qual_list  ):
		frame = 16;filtr = 10
		for i in xrange(0 , len(qual_list)-frame ):

			if len( [ x  for x in  qual_list[i:i+frame] if int(x)>=filtr ] ) ==frame:
				start = i
				break
		for i in xrange(len(qual_list)-frame, 0 ,-1  ):

			if len( [ x  for x in  qual_list[i:i+frame] if int(x)>=filtr ] ) ==frame:

				stop = i
				break
		return start,stop
	for each_f in glob.glob( path+'*.seq' ):
		seq = fasta_check( open( each_f,'rU' ) ).next()[-1].replace('\n','')   

		qual = fasta_check( open( each_f.replace('.seq','.qual'),'rU' )).next()[-1].replace('\n','')  

		start,stop = check_qual(  qual.split()   )
		trim_seq = seq[ start:stop  ]


		file_name = os.path.split(  each_f  )[-1]


		END  = open( trim_path + file_name + '.fasta'  ,'w' )


		title = '>' +file_name+'\n'
		END.write( title + trim_seq + '\n')
		QUAL = open( trim_path + file_name + '.fasta.qual'  ,'w' )
		QUAL.write(   title  )
		QUAL.write( ' '.join(  qual.split( )[ start:stop ]   ) +'\n' )
def check_path( path ):
	if os.path.exists(path):
		shutil.rmtree(path)
	os.mkdir( path )
def sweep_docu( path  , appendname):
	for file_append in appendname :
		all_data = glob.glob(path+'*.%s'%( file_append )  )

		for each_f in all_data:

			title  = each_f.split('(')[0]+'/'


			DirCreate( title   )

			shutil.move( '%s'%(  each_f  )  , title  )


def extract_html( data ):
	all_error = {}
	tag = re.escape(  '<tr bgColor="#ffffff">' )

	all_data = re.split( '''\s+(?=%s)'''%(tag),     data             )
	all_form = filter(lambda x: x.startswith('<tr bgColor="#ffffff">'), all_data)

	for each_form in all_form:
		data = re.findall( '\s+<\w+>([^<]+)</\w+' ,each_form)

		if data and data[-1] != '是':
			all_error[ data[0] ] = ''
	return all_error


def DirCreate(dir):
	if not os.path.exists(dir):
		os.mkdir(dir)

'''---------Make dir-----------'''

check_path(  outpath )
check_path(  trace_outputpath   )
check_path(  trim_path   )
check_path(  cali_path   )
check_path(   blast_path  )



Log_path = trace_outputpath+'all_log.html'
LOG = open( Log_path,'ab' )
error_record = {}
for r,c,file_list in os.walk( intputpath ):
	for each_f in file_list:
		if each_f.endswith('.zip'):
			z = zipfile.ZipFile(r+'/'+each_f, 'r')
			for f in z.namelist():
				path = os.path.join(trace_outputpath, os.path.split( f )[-1])
				if path.endswith(".html"):
					all_error = extract_html( z.read(f) )
					if all_error:
						error_record[  re.sub( '\-\d+\.zip$','',each_f  )  ] = all_error
for r,c,file_list in os.walk( os.getcwd() ):
	for each_f in file_list:
		if each_f.endswith('.zip'):
			file_tag = re.sub( '\-\d+\.zip$','',each_f  )

			z = zipfile.ZipFile(r+'/'+each_f, 'r')
			for f in z.namelist():
				file_name = os.path.split( f )[-1]
				if not f.endswith('.ab1') and not f.endswith( '.seq'  ):
					LOG.write(z.read(f))
					continue

				if file_tag not in error_record or  re.search('\(([^\)]+)',file_name).group(1) not in error_record[ file_tag ]: 

					path = os.path.join(trace_outputpath, file_name)
					if path.endswith("/"):
						DirCreate(path)

					else:

						file(path, 'wb').write(z.read(f))
				else:
					print( each_f+' \'s file named  '+file_name +' has poor quality,so it could be discard !!')

RAW = open(Log_path,'rU')
data = RAW.read()
tag = re.escape(  '<tr bgColor="#ffffff">' )
title = re.search( '''<tr bgColor="ivory">\s+(.+?)\s+</tr>''',  data,re.DOTALL).group(1)
title  = re.sub( '\s+<td noWrap>','', title)
title  = re.sub( '</td>','\t', title)
all_data = re.split( '''\s+(?=%s)'''%(tag),     data             )
all_form = filter(lambda x: x.startswith('<tr bgColor="#ffffff">'), all_data)

END = open(trace_outputpath+'all_log','w')
all_cache = {}
for each_form in all_form:

	data = re.findall( '\s+<\w+>([^<]+)</\w+' ,each_form)
	if data:
		all_cache[ '\t'.join( data ) +'\n' ] = ''

END.write( ''.join(sorted(  all_cache,key = lambda x: x.split('\t')[0]  ) ) )
all_sequencing = {}
all_ok  = {}
for each_seq  in all_cache:
	all_sequencing[  each_seq.split('\t')[1]  ] = ''

	if each_seq.split('\t' )[-1] == '是\n':
		all_ok[  each_seq.split('\t')[1]  ] = ''
ALLERROR = open(trace_outputpath+'all_error','w')
for each_seq in sorted(all_cache,key = lambda x: x.split('\t')[0]):
	if each_seq.split('\t')[1] not in all_ok:
		ALLERROR.write( each_seq )
'''-------------------------------automatic generate qual and seq file---------------------------------------------------'''		
for file_append in [ 'ab1']:
	all_data = glob.glob(trace_outputpath+'*.%s'%( file_append )  )
	need_files = [x.replace('(','\(').replace(')','\)').replace( ' ','\ ' )   for x in all_data     ]
	os.system( 'ttuner -Q -sd %s -qd %s %s'%( cali_path,cali_path, ' '.join( need_files  )  )  )

'''----------------------------sweep the document------------------------------------------'''		
check__trim( cali_path )
sweep_docu( trace_outputpath,['ab1','seq'] )
sweep_docu( cali_path,['fasta','qual'] )
sweep_docu( trim_path,['fasta','qual'] )

'''----------------------------Assembly the sequence----------------------'''

LOG = open( 'assembly_log','w'  )
root_path = os.path.abspath( trim_path )
base_path = os.path.abspath('./')
assembly_path = base_path+'/Assembly__%s'%( time.strftime('%Y-%m-%d-%H-%M-%S/',time.localtime(time.time()))  )
os.mkdir(  assembly_path )
ALL_SEQ = open( blast_path+'seq2blast.fasta','w'  )



for dirpath, subpath,filelist in os.walk(  root_path ):
	if dirpath == root_path:
		continue

	gap_name = dirpath.split('/')[-1]

	all_seq = [ x  for x in filelist if x.endswith( '.fasta' )      ]

	if len(  all_seq  ) == 0:
		LOG.write( gap_name+'\t'+'NO-reads\n' )

	else:
		os.chdir( dirpath )
		all_F = filter(lambda x: re.findall( '(F)' , x ) , all_seq)

		all_R = filter(lambda x: re.findall( '(R)' , x ) , all_seq)
		if len(all_R) < 1 :
			LOG.write( gap_name + '\t' + 'NO Reversed reads\n'  )

			combine_fasta( ALL_SEQ )
		elif len( all_F ) < 1 :
			LOG.write(  gap_name + '\t' + 'NO Foward reads\n'  )
			combine_fasta( ALL_SEQ )
		else:
			
			cache_name = '%s.fasta'%( gap_name )
			CACHE = open( cache_name,'w'  )
			QUAL = open( cache_name+'.qual','w'  )
			for each_seq in all_seq:
				CACHE.write( open(  each_seq,'rU' ).read() )
				QUAL .write(   open(  each_seq+'.qual','rU' ).read()   )
			CACHE.close()
			QUAL .close()




			assembly_log = '>%s.log'%(  gap_name  )

			os.system( ' '.join( [ 'phrap',  cache_name ,  '-new_ace', \
			                       assembly_log   
			                       ] )
			           )
			CONTIG = open( cache_name+'.contigs','rU'  )
			data = CONTIG.read()
			if not data.startswith('>'):

				LOG.write(  gap_name+'\t'+'Can not assembly\n'  )
				all_fasta( ALL_SEQ  )

			else:

				ASSEMBLY_LOG = open(  assembly_log[1:] ,'rU'  )
				contig_data_list = re.findall( '\nContig (\d+)[^\n]+\n(.+?)\n\n' , ASSEMBLY_LOG.read() ,re.DOTALL )
				contig_need = {}
				for [contig_id,contig_cont] in contig_data_list:
					if re.search( '(\S+F)', contig_cont  ) and re.search( '(\S+R)', contig_cont  ):


						LOG.write(  gap_name+'\t'+'Assembly successful!!!\n'  )

						if re.search( '((?:^|\n)C\s+.+F)',contig_cont  ):

							comp_tag = 'comp'

						else:
							comp_tag=''
						contig_need[ contig_id ] = comp_tag

						ASSEMB = fasta_check(   open( '%s.contigs'%( cache_name )  )   )

				if contig_need:
					CONTIG = open( assembly_path+ gap_name+'.contig','w')
					for t,s in ASSEMB:
						t_output= t[:-1].split('.')[0]+'-Contig\n'
						title = re.search(  'Contig(\d+)', t ).group(1)
						if title  in contig_need:
							CONTIG.write(  t_output + s  )
						ALL_SEQ.write( t_output + s + '\n'   )
					else:

						LOG.write(  gap_name+'\t'+'Can not assembly\n'  )

				else:
					
					all_fasta( ALL_SEQ  )


		os.chdir( root_path )
