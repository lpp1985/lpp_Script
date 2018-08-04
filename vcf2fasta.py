#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/7/27
from lpp import *
import multiprocessing
from optparse import OptionParser

usage = '''usage: python %prog [options] 

It can automaticly run BWA!!'''
parser = OptionParser(usage =usage ) 

parser.add_option("-r", "--ref", action="store", 
                  dest="ref",
                  type='string',
                  help="the reference seq")

parser.add_option("-v", "--vcf", action="store", 
                  dest="vcf",
                  type='string',
                  help="vcf file")

parser.add_option("-o", "--output", action="store", 
                  dest="output",
                  type='string',
                  help="output file")

(options, args) = parser.parse_args() 
REF = fasta_check( open(options.ref ,'rU' ) )
VCF = open(  options.vcf ,'rU' )
VCF.next()
VCF.next()
OUTPUT = open( options.output ,'w' )

seq = list(re.sub(  '\s+','',  REF.next()[-1]) )
print( len( seq ) )
i=0
sub_hash = {}
for line in VCF:
    i+=1
    line_l = line.split('\t')
    location = int( line_l[1]  ) 
    old = line_l[3]
    new = line_l[4].split(',')[0]
    for each_pos in xrange( location-1, location+len(old)-1   ):
        seq[  each_pos    ] = 'Y'
    sub_hash[ i ] = new
    
k=0
def subsitution(x):
    global k
    k+=1
    return  sub_hash[k]
seq = ''.join( seq )
seq_new = re.sub( 'Y+',subsitution, seq  )
seq_new = re.sub( '(\w{70})','\\1\n',seq_new  )
OUTPUT.write( '>END\n%s\n'%(  seq_new  )  )
