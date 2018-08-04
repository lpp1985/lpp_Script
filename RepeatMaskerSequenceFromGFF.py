#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/1/11
"""
import sys
from lpp import *
from optparse import OptionParser


if __name__ == '__main__':
    usage='''usage: python %prog [options]

        Transform RepeatMasker gff result'''
    parser = OptionParser(usage =usage )
    parser.add_option(
        "-i", "--Input", action="store",
        dest="fasta",
        type='string',
        help="RAW sequence")
    
    parser.add_option(
        "-g", "--gff", action="store",
        dest="gff",
        type='string',
        help="RepeatMasker gff file")    
    parser.add_option(
        "-o", "--out", action="store",
        dest="output",
        type='string',
        help="Msaked Seuqence")
    parser.add_option(
        "-s", "--seq", action="store",
        dest="seq",
        type='string',
        help="Repeat Seuqence")     
    (options, args) = parser.parse_args()
    data_hash = Ddict()
    GFF = open( options.gff,'rU'  )
    for line in GFF:
        line_l = line.split("\t")
        start = int(line_l[3])
        end = int(line_l[4])
        if start > end:
            end,start = start,end
        name = line_l[0]
        for i in xrange(start,end):
            data_hash[name][i]=""
    RESULT = open(options.output,'w')
    RAW = fasta_check(   open( options.fasta,'rU') )
    for t,s in RAW:
        name = t.split()[0][1:]
        s = re.sub("\s+","",s)
        s = list(s)
        for each_pos in data_hash[name]:
            s[each_pos]="N"
        s = "".join(s)
        s = re.sub("(\w{60})","\\1\n",s)
        RESULT.write(t+s+'\n')
