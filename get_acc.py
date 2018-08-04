#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/25
from lpp import *
RAW = open( sys.argv[1],'rU' )
RAW.next()
END = open(sys.argv[2] , 'w')
all_acc = {}
for line in RAW:
	line_l = line.split('\t')
	acc = line_l[7]
	END.write( line_l[2]+'\t'+acc+'\n' )