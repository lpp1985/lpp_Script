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
		


FASTA_ALL = fasta_check(open( intput,'rU'  ))
seq = FASTA_ALL.next()[-1]

				
