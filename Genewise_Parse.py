#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/10/20
"""

from lpp import *
from optparse import OptionParser



if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] '''
    parser = OptionParser(usage =usage )    
    parser.add_option("-s", "--Score", action="store",
                  dest="score",
                  type='int',
                  default = 60,
                  help="bitscore thread of genewise")
    parser.add_option("-k", "--Coverage", action="store",
                      dest="coverage",
                      type='float',
                      default = 0.75,
                      help="coverage thread of Protein")    



    parser.add_option("-p", "--pro", action="store",
                      dest="protein",
                      help="protein sequence!!")
    parser.add_option("-i", "--genewise", action="store",
                      dest="genewise",
                      help="Genese Result!!")
    parser.add_option("-o", "--output", action="store",
                      dest="output",
                      help="GFF Parse output!!")    
    
    (options, args) = parser.parse_args()
    score = options.score
    coverage = options.coverage
    protein = options.protein

    genewise = options.genewise
    output = options.output  
    proteinseqHash = {}
    END = open(output,'w')
    for t,s in fasta_check(open(protein,'rU')):
	if "STRG" not in t:
        	proteinseqHash[t.split()[0][1:]] = re.sub("\s+","",s) 
	else:
		name = re.search("(STRG\.\d+\.\d+)",t).group(1)
		proteinseqHash[ name ] = re.sub("\s+","",s) 
    RAW = re.split("//\nBits",open(genewise,'rU').read())
    j=0
    for e_b in RAW:
        
        data_b = e_b.split("\n//\n")
        score_b = data_b[0]
        
        aln_detail = score_b.split("\n")[1]

        aln_list = aln_detail.split()

        aln_score,alnpro,start,end = aln_list[:4]
        start = int(start)
        end = int(end)
        alnlength = end-start
        proteinlength = len(proteinseqHash[alnpro])
        
        if alnlength*1.0/proteinlength>=coverage:
            j+=1
            gff_b = data_b[1]
            for line in gff_b.split("\n"):
                if not line:
                    continue
                gff_list = line.split("\t")
                if gff_list[2] !="cds":
                    continue
                seq_name = gff_list[0]
               
                subjname,loc_append = seq_name.rsplit("__",1)
                loc_append  = int(loc_append)
                gff_list[0] = subjname
                gff_list[3] = int(gff_list[3]) + loc_append
                gff_list[4] = int(gff_list[4]) + loc_append
                gff_list[3],gff_list[4] = sorted([ gff_list[3],gff_list[4] ])
                protein_length = (gff_list[4] - gff_list[3])/3
                
                
                gff_list[3] = str(gff_list[3])
                gff_list[4] = str(gff_list[4])
                
                gff_list[-1] = "ID="+alnpro+'.%s'%(j)+';'+"Target=%s"%(alnpro+'.%s'%(j))+' %s %s'%(start,start+protein_length)
                start=start+protein_length+1
                gff_list[-2]="."
                gff_list[2] = "nucleotide_to_protein_match"
                score = 100.0*protein_length/proteinlength
                if score >100:
                    score =100
                gff_list[5]="%.2f"%(score)

                END.write("\t".join(gff_list)+'\n')
                
        
