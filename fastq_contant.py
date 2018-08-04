#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/7/12
from lpp import *
from optparse import OptionParser

usage = '''usage: python %prog [options] 

It can automaticly cut reads data!!'''
parser = OptionParser(usage =usage ) 

parser.add_option("-p", "--path", action="store", 
                  dest="path",
                  type='string',
                  help="the RAW FASTQ")

#parser.add_option("-a", "--appendix", action="store", 
                  #dest="appendix",
                  #type='string',
                  #help="the RAW FASTQ")

parser.add_option("-o", "--output", action="store", 
                  dest="output",
                  type='string',
                  help="the output File")




(options, args) = parser.parse_args() 




path = options.path

output = options.output

#appendix = options.appendix


OUTPUT = open(  output  ,'w' )

OUTPUT.write(  '#File\tReads_number\tData(MB)\n'   )

#all_files = glob.glob(   path + '*.' + appendix   )
for a1,b,c in os.walk(  path  ):
	for each_f in c:
		reads_number = 0
		RAW_raw = open( a1+'/'+each_f ,'rU'  )
		if not re.search( '(^\@)' ,  RAW_raw.next()):
			continue
		RAW = fastq_check(  open( a1+'/'+each_f ,'rU'  )   )
		data_contant = 0
		for (a,b,c,d) in RAW:
			reads_number+=1
			data_contant +=  len(b[:-1])
		OUTPUT.write( a1+'/'+each_f+'\t%s\t%.2f\n'%( reads_number,   float(  data_contant  )/1024/1024 )  )
	
	
	