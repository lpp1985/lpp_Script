#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/9
from lpp import *

from lpp import *
import multiprocessing
from optparse import OptionParser
''' sam_all_have is a hash stores all the reads name which could mapping to reference'''

sam_all_have = {}

class fastq_sam_extract(  object    ):
    def __init__( self, read1,read2  ):
        global sam_all_have
        self.READ1 = fastq_check(  open( read1,'rU' )    )
        
        self.READ2 = fastq_check(  open( read2,'rU' )    )
        self.END1 = open(  read1+'.mappable' ,'w'   )
        self.END2 = open(  read2+'.mappable' ,'w'   )
        '''t,r,qt,q represents title, readsbase,quality title, quality in the fastq file'''

        
    def __iter__(self):
        return self
    
    def next(self):
        
        
        t1,r1,qt1,q1 = self.READ1.next()
        
        t2,r2,qt2,q2 = self.READ2.next()
        name = re.search( '^@(\S+)\/',t1   ).group(1)
        if name in sam_all_have:
            self.END1.write(  ''.join( [t1,r1,qt1,q1] ) )
            self. END2.write(  ''.join( [t2,r2,qt2,q2 ])   )
        

#if __name__ =='__main__':
    
    
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'

usage='''usage: python %prog [options] 

It can stats quality and filter bad quality reads'''

parser = OptionParser(usage =usage ) 

parser.add_option("-1", "--READ1", action="store", 
                  dest="read1",
                  type='string',  
                  help="Read1 File")

parser.add_option("-2", "--READ2", action="store", 
                  dest="read2",
                  type='string',  
                  help="Read2 File")


parser.add_option("-s", "--SAM", action="store", 
                  dest="sam",
                  type='string',  
                  help="sam File")






(options, args) = parser.parse_args() 

read1 = options.read1

read2 = options.read2

sam = options.sam

if __name__ =='__main__':
    SAM_HANDLE = open( sam,'rU'   )
    all_reads_name = re.findall( '(HWUSI-[^/\s]+)()',SAM_HANDLE.read()  )
    sam_all_have = dict( all_reads_name )
    sam_class = fastq_sam_extract(  read1,read2  )
    for i in sam_class:
        pass
    