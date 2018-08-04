#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/26
"""
import sys,os
from lpp import *
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Input", action="store",
                  dest="data",

                  help=" Differential gene list by Deseq ")


parser.add_option("-l", "--Genelength", action="store",
                  dest="genelength",

                  help="All Gene Length")

parser.add_option("-g", "--Go", action="store",
                  dest="Go",

                  help="All Gene Go Annotation")

parser.add_option("-o", "--Out", action="store",
                  dest="Out",

                  help="OutputFile")

(options, args) = parser.parse_args()
data = os.path.abspath( options.data )
length = os.path.abspath( options.genelength )
go = os.path.abspath( options.Go)
GO = open( go,'rU' )
GO_Cache = open(os.path.abspath( "%s.cache"%(os.getpid()) ),'w'   )
GO_Cache.write(GO.next())
for line in GO:
    line_l = line.strip().split("\t")
    for key in line_l[1:]:
        GO_Cache.write(line_l[0]+'\t'+key+'\n')
GO_Cache.close()
end = options.Out

check_path(end)
r_script="""
library("goseq")
DEG<-read.table("%(inp)s", header = TRUE,sep="\\t")
DEG <-levels(DEG$id)
DEG.vector <- t(DEG)
ALL<-levels(read.table("%(length)s", header = TRUE,sep="\\t")[,1])
go <-read.table("%(go)s", header = TRUE,sep="\\t")

ALL.vector<-c(t(ALL))
gene.vector=as.integer(ALL.vector%in%DEG.vector)
names(gene.vector)=ALL.vector 

length_info = read.table("%(length)s", header = FALSE,sep="\\t")
gene_length.vector<-c(t(length_info[,2]))
names(gene_length.vector)<-length_info[,1]
gene_length2<-gene_length.vector[names(gene_length.vector)%in%names(gene.vector)]
pwf=nullp(gene.vector,bias.data=gene_length2,plot.fit=FALSE)
pvals <- goseq(pwf,gene2cat=go)
pvals <- pvals[ pvals$numDEInCat>1, ]
pvals$padj<-p.adjust(pvals$over_represented_pvalue, method="BH")

pvals$qvalue <-p.adjust(pvals$over_represented_pvalue, method="fdr")

pvals$gene_num<-rep(length(ALL),nrow(pvals))

pvals<-pvals[c("category","numDEInCat","numInCat","over_represented_pvalue","padj","qvalue","gene_num","term","ontology")]
colnames(pvals)[4]<-"pvalue"
pvals<-pvals[order(pvals$numDEInCat/pvals$numInCat,decreasing=T),]
#pvals<-pvals[pvals$numDEInCat/pvals$numInCat>=0.1,]
#pvals<-pvals[pvals$numInCat>=50,]
enriched_go<-pvals[pvals$qvalue<.05  ,]

write.table(enriched_go,"%(out)s/go_enrich.tsv",sep="\t",row.names=FALSE,quote=FALSE)
""".replace("%(go)s",GO_Cache.name).replace("%(inp)s",data).replace("%(length)s",length).replace("%(out)s",end)

END = open("%s.goseq.R"%(end),'w')
END.write(r_script)
END.close()

os.system("R --no-save <  %s"%(END.name))
#os.remove(GO_Cache.name)
	
