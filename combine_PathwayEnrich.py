#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/5/4
"""

from lpp import *
END = open("cache.matrix",'w')
END.write("Sample\tPathway\tRatio\tQ_value\n")
Outputhash_name = Ddict()
sample_list = {}
for f in sys.argv[1:]:
	RAW = open(f,'rU')
	RAW.next()
	for line in RAW:
		sample = os.path.split(f)[-1].split(".")[0]
		sample_list[sample] = ""
		line_l = line.strip().split("\t")
		if float(line_l[-1])<0.051:
			END.write("%s\t%s\t%s\t%s\n"%( 
			    sample,
			    line_l[1],
			    float(line_l[3])/float(line_l[2]  ),
			    line_l[-1]
			)
			        )

R_CACHE = open("Draw.R",'w')
output="enrichment_all.pdf"
R_CACHE.write("""
library(ggplot2)
countsTable <- read.delim( "%(inp)s", header=TRUE, stringsAsFactors=TRUE ) 
pathway_size = length(levels(factor(countsTable$Pathway)))
Sample_size = length(levels(factor(countsTable$Sample)))
dev.new()
pdf("%(out)s")
qplot(data = countsTable,x=Sample,y=Pathway,size=Ratio,color=Q_value)+scale_colour_gradient(low="red", high="blue")
dev.off()



"""%(
       {
           "inp" :END.name,
           "out": output
       }
      
   
   )
   
   
   )
os.system("Rscript %s"%(R_CACHE.name))
	
os.remove(R_CACHE.name)
os.remove(END.name)