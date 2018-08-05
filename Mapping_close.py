#!/usr/bin/env python
#coding:utf-8
# Author:  LPP
# Purpose: 
# Created: 2011/11/7
# Company: Chinese Human Genomic Center at Beijing
from lpp import *
from  optparse import OptionParser
usage = '''usage: python2.7 %prog [options] 

to 
'''
parser = OptionParser(usage =usage )

parser.add_option("-s", "--SCAFFOLD", action="store",
		              dest="scaff",

		              help="The scaffold to fill the gap")

parser.add_option("-o", "--OUTPUT", action="store",
                              dest="output",

                              help="The output Result")


parser.add_option("-b", "--formatdb", action="store_true", 
                  default= False,
                  dest="formatdb",
                  help="formatdb ?"
                  )

parser.add_option("-1", "--F", action="store",
                              dest="read1",

                              help="read1 fasta file built blast index")
parser.add_option("-2", "--R", action="store",
                              dest="read2",

                              help="read2 fasta file built blast index")

(options, args) = parser.parse_args()

def get_seq(   ):
	all_hash = {}
	for each_f in [ blast_read1_parse ,blast_read2_parse  ]:
		RAW = open( each_f,'rU'  )
		for line in RAW:
			line_l = line.split('\t')
			name = line_l[6][:-1]
			all_hash[name] = ''
	OUT = open( blast_read1_parse.split('.')[0]+'.allReads' ,'w')	
	for each_f in [read1,read2   ]:
		
		for a,b in fasta_check( open(each_f,'rU') ):
			name = a[1:-2]
			if name in all_hash:
				OUT.write( a+b )
	print(  OUT.name )
	return OUT.name



if __name__=='__main__':
	
	scaff = os.path.abspath( options.scaff )
	read1 = os.path.abspath( options.read1 )
	read2 = os.path.abspath( options.read2 )
	if options.formatdb:
		os.system( 'formatdb -i %s -p F'%( read1 ) )
		
		os.system( 'formatdb -i %s -p F'%( read2 ) )
	cache_docu = os.path.abspath('./Cache/')+'/' 
	if os.path.exists( cache_docu ):
		shutil.rmtree( cache_docu )
	os.mkdir( cache_docu )
	PHRAP_CACHE = open( cache_docu+'Phrap_contig.fasta' ,'w' )
	RAW = fasta_check( open(  scaff,'rU'  ) )
	os.chdir( cache_docu  )	
	for t,s in RAW:
		title = t[1:-1]
		s1 = re.sub('\s+','',s)
		i=1
		for each_Gap in re.finditer('(N+)',s1,re.I):
			[ gap_start, gap_stop ] = [ int(x) for x in each_Gap.span() ]

			[ gap_start, gap_stop ] = [  gap_start-150,gap_stop + 150    ]
			if gap_start <0:
				gap_start=0
			CACHE = open( cache_docu+'s'+title+'-gap%s'%(i)+'.splice','w'    )
			CACHE.write( '>s'+title+'-gap%s'%(i)+'\n'+s1[ gap_start: gap_stop ]+'\n' )
				
			
			blast_read1_result =  cache_docu+title+'-gap%s.read1'%(i)+'.xml'
			
			os.system( 'blastall -p blastn -i %s -d %s -m 7 -a 10 -e 1e-5 -F F >%s'%( CACHE.name,read1,blast_read1_result )  )
			
			blast_read2_result =  cache_docu+title+'-gap%s.read2'%(i)+'.xml'
			
			os.system( 'blastall -p blastn -i %s -d %s -m 7 -a 10 -e 1e-5 -F F >%s'%( CACHE.name,read2,	blast_read2_result)  )
			blast_read1_parse = blast_read1_result+'.Bparse'
			blast_parse(  open( blast_read1_result,'rU'  ),open(blast_read1_parse ,'w')  ).parse()
			
			
			
			blast_read2_parse = blast_read2_result+'.Bparse'
			blast_parse(  open( blast_read2_result,'rU'  ),open(blast_read2_parse ,'w')  ).parse()	
			all_reads_file = get_seq()
			os.system(   'phrap -vector_bound 0 -trim_start 0 -forcelevel 3 -preassemble -bandwidth 10 -repeat_stringency 0.98 -new_ace -minmatch 5 -maxmatch 20 -minscore 50 %s >see 2>&1' %( os.path.split(all_reads_file )[-1] )   )
			PHRAP_END = open( all_reads_file+'.contigs' ,'rU'   )
			PHRAP_CACHE.write( PHRAP_END.read()+'\n'  )
			i=i+1
	os.system(   'phrap -vector_bound 0 -trim_start 0 -forcelevel 3 -preassemble -bandwidth 10 -repeat_stringency 0.98 -new_ace -minmatch 5 -maxmatch 20 -minscore 50 %s >see 2>&1' %( PHRAP_CACHE.name )   )
	os.system( 'cat %s.contigs %s.singlets >Phrap_All  '%( PHRAP_CACHE.name, PHRAP_CACHE.name )   )
	os.system( 'close_gap.py -s %s -p Phrap_All -o %s -l left.fasta'%(  scaff  ,options.output ) )