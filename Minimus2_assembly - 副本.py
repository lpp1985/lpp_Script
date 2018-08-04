#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/22
from lpp import *
from optparse import OptionParser
import subprocess
import multiprocessing
usage = '''usage: python2.7 %prog -i [Input_PATH] -o [Prefix] -t [ Thread  ]'''
parser = OptionParser(usage =usage ) 
parser.add_option("-i", "--INPUT", action="store", 
	dest="input", 
	help="Input_PATH")
parser.add_option("-o", "--prefix", action="store", 
	dest="prefix", 
	default = 'Minimus_assembly',
	help="The PREFIX you want!!!!")
parser.add_option("-t", "--thread", action="store", 
	dest="thread", 
    type='int', 
    default = 50, 
	help="The Thread number  you want !!!!")
(options, args) = parser.parse_args() 
intput = options.input
PREFIX = options.prefix
thread = options.thread
if not os.path.isdir( PREFIX ):
	os.mkdir( PREFIX )
FASTA_END = open( PREFIX+'/END.fasta' ,'w' )
LOG =  open( PREFIX+'/END.log','w'  )
def minimus2(  (prefix,root,fasta_file)    ):
	yesterday = os.path.abspath(  './'  )
	name = re.search( '([^\.]+)',fasta_file ).group(1)
	docu = prefix+'/'+name+'/'
	if not os.path.isdir(  docu  ):
		os.makedirs( docu  )
		pipe = subprocess.PIPE
		os.system(     ' '.join(  ['toAmos','-s',root+'/'+fasta_file,'-o','%s%s.afg'%( docu, name ) ,'1>','%s%s.log'%( docu, name ),'2>','%s%s.error'%( docu, name )]  )     )
		number = subprocess.Popen( ['grep','-c','>',root+'/'+fasta_file ],stdout= pipe).stdout.read()[:-1]
		os.chdir(  docu  )
		os.system(    ' '.join(  ['minimus2','-D',number, name ,'2>','%s.error'%( name ),'1>','%s.log'%( name ) ]   )     )
		os.chdir( yesterday )
		


input_cache = []
for each_f in os.listdir( intput ):
	if os.path.isfile(intput+'/'+each_f):
		if each_f.endswith('.fasta' ):
			pipe = subprocess.PIPE
			number = subprocess.Popen( ['grep','-c','>',intput+'/'+each_f ],stdout= pipe).stdout.read()[:-1]
			if number =='1':
				FASTA_1 = fasta_check(  open( intput+'/'+each_f,'rU' )  )
				for t,s in FASTA_1:
					FASTA_END.write(  '>%s\n'%( each_f.replace( '.fasta','' )  ) +s    )
			else:
				input_cache.append( [PREFIX,intput,each_f] )

pool = multiprocessing.Pool(thread)
pool.map( minimus2, input_cache)
error_path = PREFIX+'/Error/'
os.mkdir( error_path  )
for root,doc,files_list in os.walk( PREFIX  ):
	for e_f in files_list:
		if e_f.endswith(  '.error'  ):
			size = os.path.getsize( root+'/'+e_f  )

			if size==0:
				FASTA = fasta_check( open( root+'/'+e_f.replace( '.error','.fasta'  ) ,'rU' ) )
				RAW = open( root+'/'+e_f.replace( '.error','.fasta'  ) ,'rU' )
				if not RAW.read():
					i=0
					RAW = fasta_check( open( intput+'/'+e_f.replace( '.error','.fasta'  ) ,'rU' ) )
					for t,s in RAW:
						name = '>'+e_f.replace( '.error',''  )+'_%s\n'%( i )
						FASTA_END.write(  name + s  )
					break
				i=0
				for t,s in FASTA:
					i+=1
					name = '>'+e_f.replace( '.error',''  )+'_%s\n'%( i )
					FASTA_END.write(  name + s  )
			else:
				shutil.copy(  root+'/'+e_f , './'   )
				os.remove(         root+'/'+e_f            )
				minimus2(  [PREFIX,intput,e_f.replace( '.error','.fasta'  )]  )
				LOG.write(  e_f.replace( '.error',''  )+'\t'+'ERROR\n'   )
				
