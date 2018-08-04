#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/1/20
"""
from lpp import *
def parse_bed(DATA):
    sample =[]
    start=[]
    width=[]
    strand=[]
    for line in DATA:
        line_l = line.strip().split("\t")
        sample.append("\""+line_l[0]+"\"")
        if int(line_l[1])<1:
            line_l[1]="1"
        start.append(line_l[1])
        width.append("%s"%(int(line_l[2])-int(line_l[1])   ))
        strand.append("\""+line_l[-1]+"\"")
    output = "GRanges(seqnames=c(%s),IRanges(start=c(%s),width=c(%s)),strand=c(%s)) "%(
        ",".join(sample),
        ",".join(start),
        ",".join(width),
        ",".join(strand),
    )  
    if not strand:
        return ""    
    return output



for a,b,c in os.walk(os.getcwd()):
    if glob.glob(a+'/*.bam'):
        SCRIPT = open(a+'/R.script','w')
        SCRIPT.write("""
library("ggbio")
library("Rsamtools")
library(GenomicAlignments)

		""")	
        all_sample = []
    for e_f in c:
        if e_f.endswith(".bam"):
            bam = a+'/'+e_f
            bed = os.popen("bamToBed -i %s"%(a+'/'+e_f))
            bed_result = parse_bed(bed)
            if not bed_result:
                continue
            sample = e_f.split(".")[0]
            all_sample.append(sample)

            SCRIPT.write("""bamfile%s<-%s \n"""%(sample,bed_result))
            SCRIPT.write("""coverage%s <- autoplot(bamfile%s,stat="coverage",fill="blue") \n"""%(sample,sample))
            SCRIPT.write("""max%s <- as.integer( max(coverage(bamfile%s)) ) \n"""%(sample,sample))
            SCRIPT.write("""align%s <- autoplot(bamfile%s,aes(fill = strand)) \n"""%(sample,sample))
            SCRIPT.write("""jpeg(file="%s.jpeg",height=1000,width=1000)\n"""%(a+'/'+sample))
            SCRIPT.write("""ff%(sample)s<-tracks("Coverage of %(sample)s"=coverage%(sample)s  ,"Alignment of %(sample)s"=align%(sample)s)\n
ff%(sample)s			
			"""%(
                               {"sample":
                                sample
                                }   
                           )
                           )
            #SCRIPT.write("""ggsave(ff%s,filename="%s.png")\n"""%(sample,a+'/'+sample))
            SCRIPT.write("dev.off()\n")
    if glob.glob(a+'/*.bam'):
        cov_cache =[]
        aln_cache=[]
        max_cache = []
        max_total="total_max<-max(c("
        for e_sample in all_sample:
            max_cache.append("max%s"%(e_sample))
        max_total+=','.join(max_cache)+'))\n'
        SCRIPT.write(max_total)
        
        for e_sample in all_sample:
            cov_cache.append(""""Coverage of %(sample)s"=coverage%(sample)s """%({"sample":e_sample}))
            aln_cache.append(""""Alignment of %(sample)s"=align%(sample)s"""%({"sample":e_sample}))
        SCRIPT.write("""jpeg(file="%s.jpeg", height=%s,width=1920)\n"""%(
            a+'/'+"total", 500*len(cov_cache)
        )
                    )
        
        SCRIPT.write(  "total<-tracks(%s,%s)+ ylim(0,1.1*total_max)\n"%(",".join(cov_cache),",".join(aln_cache)))
        SCRIPT.write("total\n")
        SCRIPT.write("dev.off()")
        os.system("nohup Rscript %s&"%(SCRIPT.name))
        #SCRIPT.write("""ggsave(total,fulename="%s.png")"""%(a+'/total'))
