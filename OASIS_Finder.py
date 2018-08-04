#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/3/13
"""

from lpp import *
import os

usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-g", "--GBK", action="store",
                  dest="genbank",

                  help="Genome Sequence in fasta format")
parser.add_option("-o", "--out", action="store",
                  dest="outputprefix",

                  help="oututprefix")



if __name__ == '__main__':
	pid = str( os.getpid() )
	(options, args) = parser.parse_args()
	gbk = options.genbank
	outputprefix = options.outputprefix
	os.system(   "OASIS  -g %s -o %s  "%( gbk,pid )    )
	os.system(   "sort -n -k 5 %s > %s"%( pid+'.gff' ,  pid+'.gff2'   ) )
	GFF = open(  pid+'.gff2' ,'rU')
	source = open(gbk).next().split()[1]
	SEQ = fasta_check(  open("%s.fasta"%(pid),'rU')  )
	all_te_seq = {}
	for t,s in SEQ:
		if "ORF" in t:
			continue
		t= t.split("|")[0]
		t_list = t.split("_")
		name = "_".join(t_list[-2:])
		all_te_seq[ name ] = s
	number =0
	RESULT_GFF = open( outputprefix+'.gff','w' )
	RESULT_XLS = open( outputprefix+'.xls','w' )
	RESULT_SEQ = open( outputprefix+'.fasta','w' )	
	RESULT_XLS.write("ID\tRef_Source\tRef_Start\tRef_End\tKind\tFunction\tRef_Frame\tSeq_Nucl_Length\tSeq_Nucleotide\tIS_Family\tIS_Group\tIS_Origin\tIS_Bitscore\tIS_Evalue\tIS_Identities\tIS_Gaps\tIS_SubjectLength\tIRL\tIRR\tLocus_Tag\n")
	for line in GFF:
		line_l= line.split("\t")
		start,end = line_l[3],line_l[4]
		SEQ_TMP = open("tmp.fa",'w')
		SEQ_TMP.write( ">tmp\n"+ all_te_seq[   "_".join([start,end])  ]  )
		SEQ_TMP.close()

		os.system(   "IsFinderAlone.py -i %s -o %s"%(  SEQ_TMP.name,  "tmp.xls" ) )
		
		if os.path.getsize( "tmp.xls" ):
			TMPRES =  open("tmp.xls")
			TMPRES.next()
			all_data = TMPRES.next()
			all_data_l = all_data.split("\t")
			family = all_data_l[5]
			group =  all_data_l[6]
			
			irl = re.search("IRL \"(\S+)\"",line).group(1)
			irr =  re.search("IRR \"(\S+)\"",line).group(1)
			if irl=="False":
				irl=""
			if irr=="False":
				irr=""
			locus_tag = re.search("locus_tag \"(\S*)\"",line).group(1) 
			
			number +=1
			line_l = line.strip().split("\t")
			line_l[0] = source
			gff_detail_l = line_l[-1].split("; ")
			is_id = "%s"%( source+"_IS"+str(number) ) 
			gff_detail_l.insert( 0,"ID \"%s\""%("ID "+is_id ) ) 
			gff_detail_l[2] = """family \"%s\""""%(  family  )
			gff_detail_l[3] = """group \"%s\""""%(  group  )
			line_l[-1] = '; '.join(   gff_detail_l )
			RESULT_GFF.write(  '\t'.join(  line_l ) +'\n')
			RESULT_XLS.write( is_id+'\t'+source+'\t'+start+'\t'+end+'\t'+all_data [:-1] +'\t'+irl+'\t'+irr+'\t'+locus_tag+'\n')
			RESULT_SEQ.write(    '>'  +  is_id  +  '\n'  +    all_te_seq[   "_".join([start,end])  ]   )
			
			
			
	
	
	os.system(  "IS_Stat.py -i %s -o %s"%(  RESULT_XLS.name,  outputprefix+'.stat') ) 