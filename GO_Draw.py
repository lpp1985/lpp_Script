#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/1/23
"""
import sys,os
from optparse import OptionParser
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--STAT", action="store",
                  dest="STAT",

                  help="GO STAT File")


parser.add_option("-o", "--Output", action="store",
                  dest="OUTPUT",

                  help="Graph Appendix")
parser.add_option("-r", "--R", action="store",
                  dest="RSCRIPT",

                  help="Rscript name")
if __name__=="__main__":
	(options, args) = parser.parse_args()
	

	commandline = """


#!/usr/local/bin/Rscript
require(ggplot2)
require(ggthemes)
library(grid)
exampleFile = "%s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE )

hei<-length(levels(countsTable$Function))
bb <-ggplot(countsTable,aes(Function,GeneNumber))
aa<-bb+facet_grid(.~Component,scales="free_x",space="free")+theme_few()+geom_bar(aes(fill=Component,position="dodge",order=Component),stat="identity")+theme(legend.position="none",axis.text.x=element_text(angle=75,hjust=1.0,size=12),strip.text.x = element_text(size=14,color="darkred",face="bold")  )+ylab("Gene Number")
#tiff("%s.tiff",width=50*hei,type="cairo")
#ggplot_build(aa)
pdf("%s.pdf",width=15)
#ggplot_build(aa)
aa
dev.off()


"""%(options.STAT,options.OUTPUT,options.OUTPUT)
	SCRIPT = open(options.RSCRIPT,'w')
	SCRIPT.write(commandline)
	SCRIPT.close()
	os.system("Rscript %s"%(options.RSCRIPT))
