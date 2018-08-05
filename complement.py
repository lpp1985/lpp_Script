#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
from lpp import *

if __name__=='__main__':
	while True:
		raw_data = raw_input('Please Input Your Sequence!!!\n')
		print( complement(raw_data) )
