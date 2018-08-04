#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2014/1/20

import os,sys
nucleo = sys.argv[1]
pro = sys.argv[2]
PRO = open( pro,'w' )
data = os.popen('transeq -sequence %s -outseq stdout'%( nucleo )).read()
data = data.replace('_1\n','\n')
PRO.write(data)
