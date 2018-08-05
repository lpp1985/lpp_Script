#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/1/23
"""
import sys,os

commandline = """

require(ggplot2)
library(grid)
countsTable <- read.delim( "%s", header=TRUE, stringsAsFactors=TRUE )
countable3 <- countsTable[countsTable$Q_value<0.05,]
#countable3<-countable2[order(countable2$Q_value),]
aa<-ggplot(countable3)+geom_bar(aes(x=Name,y=Diff, fill=Q_value),stat="identity")+coord_flip()+ylab("Differential Gene Number")+theme(axis.text.y=element_text(color="darkred",face="bold"))
pdf("%s.pdf",width=10,height=3)
aa <- ggplot_gtable(ggplot_build(aa))
grid.draw(aa)
dev.off()


"""%(sys.argv[1],sys.argv[2])
SCRIPT = open(sys.argv[3],'w')
SCRIPT.write(commandline)
SCRIPT.close()
os.system("Rscript %s"%(sys.argv[3]))
