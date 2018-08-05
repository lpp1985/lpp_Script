#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from lpp import *
import subprocess
import os,sys
from os.path import abspath
from  termcolor import colored
from optparse import OptionParser
from ConfigParser import ConfigParser
import pandas as pd

general_config = ConfigParser()
path =os.path.split(__file__)[0]
general_config.read(
    path+"/general.ini"
)

def Config_Parse():
	config_hash = Ddict()
	for section in general_config.sections():
		for key in general_config.options(section):

			config_hash[section][key] = general_config.get(section,key)

	return config_hash



def Nul_or_Protein( seq ):
	seq = re.sub("\s+","",seq.lower())
	if len(set(seq) )<7 and set(seq) & set([ 'a','t','c','g'     ]):
		blast_type="blastx"
	else:
		blast_type = "blastp"
	return blast_type

def RunDiamond( fasta,evalue,blasttype,dbname,output  ):
	if not output.endswith(".tsv") and not output.endswith(".xls"):
		output=output+".tsv"
	command = "diamond_align.py  -i %(fasta)s  -n %(dbname)s -o %(output)s  -a  %(dbname)s   -e %(evalue)s  -t %(blasttype)s"%(
	    {
	        "fasta":fasta,
	        "output":output,
	        "dbname":dbname,
	        "evalue":evalue,
	        "blasttype":blasttype
	        
	    
	    
	    }
	    
	)
        print(command)
	command_list = command.split()
	diamond_process = subprocess.Popen( command_list,stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
	stdout,stderr = diamond_process.communicate()
	#return stderr
