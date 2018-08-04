#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/21
"""
import sys,shlex,os,subprocess

sys.path.append(os.path.split(__file__)[0]+'/../Lib/')
from lpp import *
from optparse import OptionParser
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Input", action="store",
                  dest="Input",

                  help="Input Fasta Sequence")
parser.add_option("-o", "--Output", action="store",
                  dest="output",

                  help="Output File")

parser.add_option("-d", "--Database", action="store",
                  dest="Database",

                  help="Database Location")
parser.add_option("-e", "--Evalue", action="store",
                  dest="Evalue",

                  help="Evalue !")
if __name__=="__main__":
    (options, args) = parser.parse_args()
    Input = options.Input
    OUTPUT = options.output
    Database = options.Database
    Path = os.path.split(OUTPUT)[0]
    if not Path:
	Path="./"
    if not os.path.exists(Path):
        os.makedirs(Path)
    E_value = options.Evalue

    RAW = fasta_check(open(Input,'rU'))
    #all_seq = {}
    #for t,s in RAW:
    #   all_seq[t[1:-1]] = s.strip()+'\n'

    commandline = """    ghostz aln  -i  %s -b 1  -a 64  -q d  -d %s   -o /dev/stdout | awk 'NR > 5 {print $0}' |awk -F'[\t]'   '$11<%s'|sed    "s/^finished$//g"   | sed  '/^$/d'   >%s 2>/dev/null"""%( 
        Input,
        
        Database,
	E_value,
        OUTPUT 
    )
    os.system(commandline)
    
