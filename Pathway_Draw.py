#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/1/23
"""
import sys,os
import pandas as pd
from optparse import OptionParser

usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Input", action="store",
                  dest="input",

                  help="input file")
parser.add_option("-o", "--output", action="store",
                  dest="Output",

                  help="Output Prefix")
parser.add_option("-r", "--R", action="store",
                  dest="R",

                  help="R File")


if __name__ == '__main__':
	(options, args) = parser.parse_args()

	commandline = """
#!/usr/local/bin/Rscript
require(ggplot2)
require(ggthemes)
library(grid)
exampleFile = "%s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE )
countsTable <- countsTable[ countsTable$Category!='' ,   ] 
aa<-ggplot(countsTable)+geom_bar(aes(x=Category, fill=Category),show_guide =FALSE)+coord_flip()+ylab("Gene Number")+theme_few()+theme(axis.text.y=element_text(size=15,color="darkred",face="bold"))
#tiff("%s.tiff", width=1024, height=512,type="cairo")
#ggplot_build(aa)

pdf("%s.pdf",width=15)
#ggplot_build(aa)
aa
dev.off()



	
	
	"""%(options.input,options.Output,options.Output)
	SCRIPT = open(options.R,'w')
	SCRIPT.write(commandline)
	SCRIPT.close()
	r_path = os.path.split( os.path.abspath(options.R )  )[0]
	os.system("Rscript %s &&rm %s/Rplot*.pdf"%(options.R,r_path))
