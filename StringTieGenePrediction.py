#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/10/19
"""
from lpp import *
import os
from optparse import OptionParser
def check_path( path ):
	if not os.path.exists(path):

		os.makedirs( path )
	return os.path.abspath(path)+'/'


if __name__=='__main__':
	usage = '''usage: python2.7 %prog [options] Kmer




    Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
	parser = OptionParser(usage =usage )





	parser.add_option("-g", "--gtf", action="store",
                      dest="gtf",
                      help="string gtf Sequence!!")

	parser.add_option("-a", "--assembly", action="store",
                      dest="assembly",
                      help="Assemblied Genome!!")


	parser.add_option("-t", "--pep", action="store",
                      dest="Transdecoder",
                      help="TransDecoder protein seqeunce")    

	parser.add_option("-o", "--out", action="store",
                      dest="output",
                      default = 'genewise.out',
                      help="The output file  you want!!")


	(options, args) = parser.parse_args()

	gtf = options.gtf
	protein = options.Transdecoder
	assembly = options.assembly
	output = os.path.abspath(options.output)

	assemblyseqHash = {}
	for t,s in fasta_check(open(assembly,'rU')):
		t  = t.split()[0][1:]
		s = re.sub("\s+",'',s)
		assemblyseqHash[t]=s

	est_hash = {}
	GFF = open(gtf, 'rU')
	#GTF infomation parse
	for line in GFF:
		if "#" in line:
			continue
		line_l = line.split("\t")
		if line_l[2] == "transcript":
			scaf_name = line_l[0]
			transcript_id = re.search("""transcript_id \"(STRG[^\"]+)\"""", line).group(1)

			scaffoldStart,scaffoldEND = int(line_l[3]), int(line_l[4])
			direction =  line_l[6]
			if scaffoldStart<10000:
				scaffoldStart = 0
			else:
				scaffoldStart =scaffoldStart -10000
			if scaffoldEND + 10000 > len(assemblyseqHash[scaf_name]):
				scaffoldEND =  len(assemblyseqHash[scaf_name])
			est_hash[transcript_id] = [scaf_name, scaffoldStart,scaffoldEND, direction ]



	COMMAND = open("genewise.cmd",'w')

	cache_path = check_path("CACHE/")
	i=0
	for t,s in fasta_check(open(protein,'rU')):

		if "complete" in t:
			t = t.strip()


			i += 1
			proteinid = re.search("(STRG[^\:]+)", t).group(1)
			PRO= open(cache_path+'%s.pep'%(i),'w')
			PRO.write('>'+proteinid+'\n'+s+'\n')            
			NUC = open(cache_path+'%s.nuc'%(i),'w')
			est_info = est_hash[ proteinid]
			scaf_name, scaffoldStart,scaffoldEND,direction =  est_info
			NUC.write( ">"+scaf_name+"__%s\n"%(scaffoldStart) )
			NUC.write( assemblyseqHash[scaf_name][scaffoldStart:scaffoldEND] + '\n')

			if direction == "+" and  t[-2] == "+" or   direction == "-" and  t[-2]  == "-":
				commandline = """genewise %s  %s   -gff -silent -quiet  -sum    """%(PRO.name,NUC.name)
			else:
				commandline = """genewise %s  %s   -gff -silent -quiet  -sum  -trev   """%(PRO.name,NUC.name)


			COMMAND.write(commandline + '\n')
	COMMAND.close()
	os.system("Genewise.nf --command %s --out %s" % ( COMMAND.name , output) )

