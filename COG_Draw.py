#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/1/23
"""
import os,sys
from os.path import abspath
from optparse import OptionParser
from lpp import *

usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--CAT", action="store",
                  dest="CAT",

                  help="COG CAT File")


parser.add_option("-o", "--Output", action="store",
                  dest="OUTPUT",

                  help="Graph Appendix")
parser.add_option("-r", "--R", action="store",
                  dest="RSCRIPT",

                  help="Rscript name")

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    commandline = """
#!/usr/local/bin/Rscript
require(ggplot2)
require(ggthemes)
library(grid)
exampleFile = "%s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE )
countsTable = countsTable[countsTable$COG_FunCat!='',]
aa<-ggplot(countsTable)+geom_bar(aes(x=COG_Category.Annotation, fill=COG_Category.Annotation),show_guide =FALSE)+coord_flip()+ylab("Gene Number")+theme_few()
pdf("%s.pdf",width=15  )
#ggplot_build(aa)
aa
#tiff("%s.tiff", width=1024, height=512,type="cairo")
#ggplot_build(aa)
dev.off()


"""%( options.CAT,options.OUTPUT, options.OUTPUT )
    SCRIPT = open(options.RSCRIPT,'w')
    SCRIPT.write(commandline)
    SCRIPT.close()
    os.system("Rscript %s"%(options.RSCRIPT))
