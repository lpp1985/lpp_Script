#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/11/5
"""

from lpp import *
import itertools
all_pep = glob.glob(sys.argv[1]+"/*.pep")
PSUDEO = open()
pep_length = {}
for e_pep in all_pep:
    os.system(  "makeblastdb  -in %s -title ref  -parse_seqids  -out %s -dbtype nucl "%(e_pep,e_pep)  )
    PEP = fasta_check(open(e_pep,'rU'))
    for t,s in PEP:
        s = t.strip()[:-1]
        t = t.strip()[1:].split()[0]
        pep_length[t] = len(s)

blast_data = Ddict()
psuedo = Ddict()
if __name__ == '__main__':
    ORTHOLOG = open("OUT/1.Orthologs_Cluster.txt",'rU')
    PSUEDO = open("PsuedoGene.list",'w')
    PSUEDO.write("Query_name\tQuery_Kind\tSubject_name\tSubj_Kind\tSubject_percent\tScore	E_value	Positive	Identity	Query_length	Query_start	Query_end	Query_percent	Subject_length	Subject_start	Subject_end	Subject_description\n")
    ORTHOLOG.next()
    for line in ORTHOLOG:
        line_l = line.strip().split("\t")[1:]
        most_length_data = sorted(line_l,key= lambda x: pep_length[x])[-1]
        longest = pep_length[  most_length_data ]
        for e_pep in line_l:
            if pep_length[e_pep]*100.0/longest<0.75:
                psuedo[e_pep][most_length_data] = ""
            
    
    all_comb = itertools.permutations(all_pep,2)
    for query_name,subject_name in all_comb:
        q_name = os.path.split(query_name)[-1].split(".")[0]
        s_name = os.path.split(subject_name)[-1].split(".")[0]
        name = "%s_%s"%( q_name,s_name  )
        os.system(  "blastp   -query %s -db %s  -max_target_seqs  1 -evalue 1e-10  -outfmt 5  -num_threads  64 -out %s.xml 2>/dev/null "%(query_name,subject_name,name)  )
        os.system(  "blast_pasrse.py %s.xml"%(name)  )
        RAW = open("%s.Bparse"%(name),'rU')
        END.write("Query_name\tQuery_Kind\tSubject_name\tSubj_Kind\tSubject_percent\tScore	E_value	Positive	Identity	Query_length	Query_start	Query_end	Query_percent	Subject_length	Subject_start	Subject_end	Subject_description\n")
        END = open("%s.blast"%(name),'w')
        RAW.next()
        for line in RAW:
            line_l = line.strip().split("\t")
            q_kind = "CDS"
            s_kind = "CDS"
            query= line_l[2]
            query_length = line_l[3]
            query_start = line_l[13]
            query_end = line_l[14]
            sub_start,sub_end = line_l[15],line_l[16],
            subj = line_l[5]
            subj_length = line_l[8]
            align_perc = int(line_l[-4])*100.0/int(subj_length)
            if query in psuedo:
                q_kind = "pseudo_gene"
            if subject in psuedo:
                s_kind = "pseudo_gene"
            postive = "(%s/%s) %.2f"%(line_l[-5],subj_length  , int(line_l[-5])*100.0/int(subj_length))
            identity = "(%s/%s) %.2f"%(line_l[-6],subj_length  , int(line_l[-6])*100.0/int(subj_length))
            query_coverage = "(%s/%s) %.2f"%(line_l[-4],query_length  , int(line_l[-4])*100.0/int(query_length))
            subject_coverage =  "(%s/%s) %.2f"%(line_l[-4],query_length  , int(line_l[-4])*100.0/int(subj_length))
            postive_detail = "%s/%s %.2f"%( line_l[21],subj_length,align_perc    )
            bitscore= line_l[10]
            e_value = line_l[12]
            data_l = [query,q_kind,subj,s_kind,subject_coverage,bitscore,e_value,postive_detail,identity,query_length,query_start,query_end,query_coverage,subj_length,sub_start,sub_end]
            output = '\t'.join( data_l )
            END.write(output+'\n')
            if query in psuedo:
                if subject in psuedo[query]:
                    PSUDEO.write(output+'\n')
            blast_data[query][subj][align_perc]= output 
            