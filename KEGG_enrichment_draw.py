#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/27
"""
import os
from lpp import *
ALL_Data = open(sys.argv[1])
TMP = open("%s.tmp"%(os.getpid()),'w')
TMP.write(ALL_Data.next()[:-1]+'\tSituation\n')

path = sys.argv[-1]
    
for e_f in sys.argv[1:-1]:

    sample_name = os.path.dirname(e_f).split('/')[-1]
    RAW = open(e_f,'rU')
    RAW.next()
    for line in RAW:
        TMP.write(line[:-1]+'\t'+sample_name+'\n')
    
    
    

R = open("%s.R"%(os.getpid()),'w')
r_script = """
library(ggplot2)
require(ggthemes)

go_data <- read.delim( "%(input_data)s", header=TRUE, stringsAsFactors=TRUE ) 
height = length( levels( go_data$Name  )  )/2
if (height<10) height<-10
go_data$EnrichFactor <- go_data$Diff_In / go_data$All_In
pdf("%(path)s/KEGGEnrich.pdf",width=15,height= height)

p <- qplot(Situation, Name, data=go_data, size=EnrichFactor,color=Q_value)
p +scale_colour_gradient(low="red", high="blue")+theme_few()

dev.off()



"""%(
       {
           "input_data":TMP.name,
           "path":path
      
       }
   )
R.write(r_script)
os.system("Rscript %s"%(R.name))
#os.remove(R.name)
#os.remove(TMP.name)
