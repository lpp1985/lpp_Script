#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2012/10/24

from lpp import *

from optparse import OptionParser
if __name__=='__main__':
	usage = '''usage: python2.7 %prog [options] Kmer




	Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
	parser = OptionParser(usage =usage )



	parser.add_option("-i", "--INPUT", action="store",
	                  dest="input_file",

	                  help="Input FASTA Sequence")
	
	parser.add_option("-p", "--PERFIX", action="store",
	                  dest="output_prefix",

	                  help="Out prefix")
	
	parser.add_option("-t", "--THREDSHOLD", action="store",
	                  dest="threshold",
	                  type="int",
	                  help="length to tolerate gap")
	
	parser.add_option("-b", "--BASE", action="store",
		                  dest="base_score",
		                  type="int",
		                  help="the base score to each site")	
	parser.add_option("-g", "--GAP", action="store",
		                          dest="gap_file",

		                          help="file to store gap")
	
	parser.add_option("-n", "--NAME", action="store",
		                          dest="name",

		                          help="name of end sequence")	
	
	parser.add_option("-j", "--JUMP", action="store",
	                  dest="jump_threshold",
	                  type="int",
	                  help="threshold of region between low coverage region")	
	(options, args) = parser.parse_args()
	input_file = options.input_file
	output_prefix = options.output_prefix
	threshold = options.threshold
	base_score = options.base_score
	gap_file = options.gap_file
	jump_threshold = options.jump_threshold
	name = options.name
	
	
	
	FASTA = fasta_check(  open( input_file,'rU' )   )
	GAP_FILE =  open( gap_file,'rU' )   
	
	
	gap_hash = {}
	low_quality_hash = {}
	
	cache = int(GAP_FILE.next().split()[0] )
	start = cache
	for line in GAP_FILE:
		end = int(line.split()[0])
		if (end - cache) > jump_threshold:
			if (cache-start) >threshold:
				for i in xrange( start,cache+1 ):
					gap_hash[i] = ''
				print( cache,start )
					
			else:
				for i in xrange( start,cache+1 ):
					low_quality_hash[i] = ''
			start = end
		cache = end
		
	for each_site in low_quality_hash:
		if each_site in gap_hash:
			print( each_site )
	seq = re.sub( '\s+','',FASTA.next()[-1] )
	quality_matrix  = ''
	i=0
	end_seq = ''
	for each_site in seq:
		
		i+=1
		
		if i in low_quality_hash:
			quality_matrix+= chr(65)
			end_seq+=each_site
		elif i in gap_hash:
			end_seq+='@'
			quality_matrix+='@'
		else:
			end_seq+=each_site
			quality_matrix+= chr(base_score+64)
			
END_SEQ = open(  output_prefix+'.fasta','w'  )

END_QUAL = open(  output_prefix+'.qual','w'  )


seq_list = re.split( '\@+',   end_seq   )

qual_list =  re.split( '\@+', quality_matrix   ) 


print( len(  qual_list) )
print( len(  seq_list) )
i=0
for each_sub in seq_list:
	i+=1
	END_SEQ.write( '>%s_%s\n'%( name,i  )  )
	
	sequence = re.sub( '(\w{60})','\\1\n', each_sub )
	
	END_SEQ.write( sequence+'\n' )
	
	END_QUAL.write( '>%s_%s\n'%( name,i  )  )
	

	
	quality_for_number_list  = [  str(ord( x )-64) for x in qual_list[i-1]   ]
	
	j=0
	for each_qual in quality_for_number_list:
		j+=1
		if j%60 ==1:
			END_QUAL.write(  each_qual )
		elif j%60==0:
			END_QUAL.write(  ' '+each_qual+'\n' )
		else:
			END_QUAL.write(  ' '+each_qual )
	END_QUAL.write('\n')
	
	