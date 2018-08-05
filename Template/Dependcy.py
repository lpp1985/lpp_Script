#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from lpp import *
from ConfigParser import ConfigParser


general_config = ConfigParser()
path =os.path.split(__file__)[0]
general_config.read(
    path+"/config.ini"
)

def Config_Parse():
	config_hash = Ddict()
	for section in general_config.sections():
		for key in general_config.options(section):

			config_hash[section][key] = general_config.get(section,key)

	return config_hash


class fasta_check(object):
	def __init__(self,file_handle):
		assert isinstance( file_handle,file ),'The Input paramater must be a File Handle'
		self.file=iter(file_handle)
		for line in self.file:
			if line[0]=='>':
				self.define=line
				break
	def __iter__(self):
		return self
	def next(self):
		if not self.define:
			raise StopIteration

		name=self.define
		self.define=''        
		s=[]
		for line in self.file:
			if line[0]!='>':
				s.append(line)
			else:
				self.define=line
				break
		s=''.join(s)
		return (name,s)
def check_path(path):
	path = os.path.abspath(path)
	if not os.path.exists(path):
		os.makedirs( path )
	return path+'/'
