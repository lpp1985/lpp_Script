#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/19
from lpp import *
import string
data_f = sys.argv[1:]

for e_f in data_f:
	RAW = blast_parse(open(e_f),open(e_f.rsplit('.',1)[0]+'.Bparse','w'))
	RAW.parse()
	
