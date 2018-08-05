#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Transfer Cutadpter result to fastq
# Created: 2011/6/2
from lpp import *
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly transfer cutadpt result to fastq which fetch to miRexpress to analysis'''
parser = OptionParser(usage =usage ) 

parser.add_option("-i", "--INPUT", action="store", 
                  dest="input",
                  type='string',  
                  help="the out of cutadpt")

parser.add_option("-u", "--UN", action="store", 
                  dest="un",
                  type='string',  
                  help="the untrim outupt of cutadpt")

parser.add_option("-o", "--OUT", action="store", 
                  dest="out",
                  type='string',  
                  help="the untrim outupt of cutadpt")


(options, args) = parser.parse_args() 
def fastq( RAW  ):
	for line in RAW:
		if line.startswith('@'):
			a = line
			i=0
		else:
			i+=1
			if i==1:
				b = line
			if i==2:
				c = line 
			if i==3:
				d = line
				yield a,b,c,d


END = open( options.out,'w'   )

RAW_input = fastq( open( options.input , 'rU'  )    ) 

RAW_notrim = fastq ( open( options.un , 'rU'  )     )
RAW = open( options.input , 'rU'  ) 
for a,b,c,d  in RAW_input :

	if len( re.search( '(\S+)',b  ).group(1)   )<17:
		continue
	else:
		c = c[0] +a[1:]
		END.write( a+b+c+d  )
		
for a,b,c,d in RAW_input:
	if len( re.search( '(\S+)',b  ).group(1)   )<17:
		continue
	else:
		c = c[0] +a[1:]
		END.write( a+b+c+d  )

#'''The last 15 base is droped for the reason of reliability'''	
		
for a,b,c,d in RAW_notrim:
	b = b[:-16]+'\n'
	d = d[:-16]+'\n'
	c = c[0] +a[1:]
	END.write( a+b+c+d  )
	