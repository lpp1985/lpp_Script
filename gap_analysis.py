#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/28
from lpp import *
import shutil
from optparse import OptionParser 
usage = '''usage: %prog -i Sequence_with_gap -o Output'''
parser = OptionParser(usage =usage )
parser.add_option("-i", "--SEQUENCE", action="store", 
                  dest="intput", 
                  help="Fasta file with gap represented by N")
parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output", 
                  help="Output Document")
(options, args) = parser.parse_args() 
seq_file   = options.intput
output = options.output
def check_path( path ):
	if os.path.exists(path):
		shutil.rmtree(path)
	os.makedirs( path )
RAW = fasta_check(  open(seq_file,'ru')  )

critique = 600
def check( length,scaf_length  ):
	if length<0:
		length =0
	elif length>scaf_length:
		length= scaf_length
	return length
for t,s in RAW:
	title =re.search( '>(\S+)',t ).group(1)
	path = output+'%s/'%( title )
	cand_seq_path = path+'CAN/CAN_FASTA/'
	no_seq_path = path+'NO/NO_FASTA/'
	check_path(path)
	check_path(cand_seq_path)
	check_path(no_seq_path)
	i=0
	LONG = open( path+title+'.2LONG','w'  )
	LONG.write(  'SCARF\tGAP_START\tGAP_END\tLENGTH\n')
	CANDIT = open( path+title+'.CANDIT','w'  )
	CANDIT.write(  'SCARF\tGAP_START\tGAP_END\tLENGTH\n')
	CANDIT_LOG = open(  path+title+'.PCR','w'  )
	
	s = re.sub('\s+','',s)
	length = len(  s  )
	all_s_location = re.finditer( '(N+)' ,s )
	for s_location in all_s_location:
		i = i + 1
		gap_start, gap_end = s_location.span()
		
		len_gap = gap_end - gap_start
		if len_gap >critique-200:
			LONG.write( t[:-1]+'\t'+str( gap_start )+'\t'+str( gap_end )+'\t'+str(len_gap)+'\n' )
		else:
			pen_start = (critique-len_gap)/2
			CANDIT.write  ( '%s\t%s\t%s\t%s\n'%( t[:-1],gap_start , gap_end,  len_gap  ) )
			log_start = check(  gap_start - pen_start ,  length  )
			log_end = check(  gap_end + pen_start ,  length  )
			gap_location = (gap_start - log_start+1 ,gap_end- log_start+1)
			CANDIT_LOG.write( '%s\t%s\t%s\n'%(title,log_start,log_end ) )
			CAN_GAP_FASTA = open( cand_seq_path+'gap%s.fasta'%( i ),'w' )
			cand_seq = re.sub('(\w{60})','\\1\n', s[ log_start : log_end   ] )
			CAN_GAP_FASTA.write(  '>'+t[:-1]+'|gap%s\n'%(i)+cand_seq  )
			CAN_DETAIL = open( path+'CAN/'+'gap%s.log'%( i ),'w' )
			CAN_DETAIL.write(  'No.\tgap_from\tgap_end\tgap_length\n' )
			CAN_DETAIL.write('gap%s\t%s\t%s\t%s\n'%( i,gap_location[0],gap_location[1],len_gap   ))