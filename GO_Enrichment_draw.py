#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/27
"""
import os
from lpp import *
ALL_Data = open(sys.argv[1],'rU')
UP = open( sys.argv[2],'rU')
DOWN = open(sys.argv[3],'rU')
TMP = open("tmp",'w')
TMP.write(ALL_Data.next()[:-1]+'\tSituation\n')
sample_name = os.path.split(sys.argv[1])[-1].split('.')[0]
for line in ALL_Data:
    TMP.write(line[:-1]+'\t'+sample_name+'\n')
    


sample_name = os.path.split(sys.argv[2])[-1].split('.')[0]
UP.next()
for line in UP:
    TMP.write(line[:-1]+'\t'+sample_name+'\n')
    
    
    
sample_name = os.path.split(sys.argv[3])[-1].split('.')[0]
DOWN.next()
for line in DOWN:
    TMP.write(line[:-1]+'\t'+sample_name+'\n')
    
R = open("GO_EnrichmentDraw.R",'w')
r_script = """
library(ggplot2)
require(ggthemes)
library(stringr)


go_data <- read.delim( "%(input_data)s", header=TRUE, stringsAsFactors=TRUE ) 
go_data$EnrichFactor = go_data$numDEInCat/go_data$numInCat
height = length(levels(go_data$category))
width = 5 *length(levels(go_data$Situation))
pdf("GOEnrich.pdf",width=width,height=0.2*height )
countsTable$Term <- str_wrap(countsTable$Term, width = 40)
p <- qplot(Situation, term, data=go_data, size=EnrichFactor,color=qvalue)
p + scale_size("EnrichFactor")+scale_color_gradient(low="red", high="blue")+theme_few()+facet_grid(.~ontology,scales="free_x",space="free")+theme(axis.text.x=element_text(angle=75,hjust=1.0),strip.text.x = element_text(size=14,color="darkred",face="bold"))

dev.off()



"""%(
       {
           "input_data":TMP.name
      
       }
   )
R.write(r_script)
os.system("Rscript %s"%(R.name))
