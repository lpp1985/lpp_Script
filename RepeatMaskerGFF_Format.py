#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/11
"""
import sys
from optparse import OptionParser


if __name__ == '__main__':
    usage='''usage: python %prog [options]

        Transform RepeatMasker gff result'''
    parser = OptionParser(usage =usage )
    parser.add_option(
        "-i", "--Input", action="store",
        dest="out",
        type='string',
        help="RepeatMasker out file")
    
    parser.add_option(
        "-g", "--gff", action="store",
        dest="gff",
        type='string',
        help="RepeatMasker gff file")    
    parser.add_option(
        "-o", "--out", action="store",
        dest="output",
        type='string',
        help="Output Gff File")     
    (options, args) = parser.parse_args()
    data_hash = {}
    GFF = open( options.gff,'rU'  )
    REPOUT = open( options.out,'rU'  )
    RESULT = open( options.output,'w')
    REPOUT.next()
    REPOUT.next()
    REPOUT.next()
    for line in REPOUT:
        line = line.strip()
        line = line.lstrip()
        line_l = line.split()
        data_hash[ "\t".join(  [line_l[4] ,line_l[5],line_l[6]] ) ] = line_l[10]
    for line in GFF:
        if line.startswith("#"):
            RESULT.write(line)
            continue
        line_l = line.strip().split("\t")
        line_l[-1]= data_hash[ "\t".join(  [line_l[0] ,line_l[3],line_l[4]] ) ] 
        RESULT.write(  "\t".join(line_l)+'\n'   )
            
