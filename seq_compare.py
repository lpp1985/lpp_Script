#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 321
# Created: 2011/11/7
from lpp import *
import glob,re,sys
from optparse import OptionParser 
usage = ''' To do Desq of Reads mapping matrix !
Usage python %prog     [Options]''' 
parser = OptionParser(usage =usage ) 
parser.add_option("-r", "--RAW", action="store", 
                  dest="raw_data",
                  help="raw_data")

parser.add_option("-n", "--new", action="store", 
                  dest="new", 

                  help="New data")

parser.add_option("-o", "--Output", action="store", 
                  dest="output", 
                  default = 'Compare.comp',
                  help="Compare end data")

if __name__=='__main__':
    (options, args) = parser.parse_args()
    already = {}
    raw_data_hash = {}
    RAW_UNIQUE = open(  options.raw_data+'.unique'  , 'w' )
    NEW_UNIQUE = open(  options.new+'.unique'  , 'w' )
    COMPARE = open(  options.output,'w'  )
    COMPARE.write('%s\t%s\n'%(  options.raw_data, options.new  )  )
    RAW = fasta_check( open( options.raw_data,'rU' )  )
    NEW = fasta_check( open( options.new,'rU' )  )
    for t,s in RAW:
        s = re.sub('\s+','',s  )
        raw_data_hash[s] = t[1:-1]
    for t1,s1 in NEW:
        s1 = re.sub( '\s+','',s1 )
        
        if s1 in raw_data_hash:
            already[s1] = ''

            COMPARE.write( raw_data_hash[s1]+'\t'+t1[1:-1]+'\n' )
        else:
            NEW_UNIQUE.write( t1+s1+'\n' )
            COMPARE.write( '-\t'+t1[1:-1]+'\n' )
            
for key1 in raw_data_hash:
    if key1 not in already:
        COMPARE.write( raw_data_hash[key1]+'\t-\n')
        RAW_UNIQUE.write(  '>'+raw_data_hash[key1]+'\n'+key1+'\n')
        