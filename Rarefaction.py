#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/2/3
"""

from lpp import *
import collections,sys,HTSeq
def Output(FILE,count):
    END = open(FILE,'w')
    for gene_id in counts:
        END.write('%s\t%s\n'%(gene_id, 10**11*float(counts[ gene_id ]  )  /seq_length[gene_id]/count )   ) 
        
        



if sys.argv[1].endswith(".sam"):
    bam_file = sys.argv[1].replace(".sam",".bam")
    os.system("samtools -bS %s >%s"%(  sys.argv[1],bam_file ))
else:
    bam_file = sys.argv[1]
bam_file = os.path.abspath(bam_file)
counts = collections.Counter( )
all_Data = os.popen("bamtools stats -in %s "%(bam_file))
total_reads = int(re.search("Total reads\:\s+(\d+)",all_Data.read()).group(1))/2
#total_reads = 33558914/2



fasta_file = sys.argv[2]
fasta = fasta_check(open(fasta_file,'rU'))
seq_length = {}
for t,s in fasta:
    seq_length[t.split()[0][1:]] = len(re.sub("\s+", "", s))


splice  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9]
sub_reads = [x * total_reads for x in splice ]


splice = iter(splice)

out_hash = {}
for i in splice:
    out_hash[total_reads*i] = re.sub("\.bam$","-"+str(int(i*100))+"%.count",bam_file) 
    
out_name = {}
has = {}

min_data = min(out_hash)

almnt_file = HTSeq.BAM_Reader( bam_file )

for almnt in almnt_file:
    if  almnt.aligned and almnt.proper_pair:

        gene_name = almnt.iv
        gene_name=gene_name.chrom

        counts[ gene_name ] += 1  
        
    has[almnt.read.name] = ""
    if len(has)>min_data and out_hash:
        Output(out_hash[min_data],min_data*2)
        del out_hash[min_data]
        if out_hash:
            min_data = min(out_hash)    


    
Output(re.sub("\.bam$","-100%.count",bam_file),total_reads*2  )
count_path = os.path.dirname(bam_file)+"/*.count"
os.system("matrix_combine.py %s"%(count_path))