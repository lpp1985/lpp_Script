#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/4/29
"""
from lpp import *
import os
#usage python2.7 expression_matrix_build.py     [ count.py's end__1_append_file_name  ]  [   matrix_end   ]
import glob,re,sys
from optparse import OptionParser 
usage = ''' To do Desq of Reads mapping matrix !
Usage python %prog     [Options]''' 
parser = OptionParser(usage =usage ) 
parser.add_option("-o", "--OUTPUTPATH", action="store", 
                  dest="outputpath",
                  default = 'Static_output', 
                  help="OUTPUTPATH")

parser.add_option("-a", "--append", action="store", 
                  dest="append", 
                  default = 'matrix',
                  help="Matrix of Reads number")



parser.add_option("-i", "--INPUT", action="store", 
                  dest="input_path", 

                  help="input path of expression")


parser.add_option("-g", "--GOMAPPING",  
                  dest="gomapping", default=True,  
                  help="A file to store gene and GO mapping!") 



(options, args) = parser.parse_args()

outputpath   = options.outputpath


input_path   = options.input_path
 
append = options.append

gomapping = options.gomapping
def sampleNameTrans( sample_name ):
	if re.search( '(^\d+)',sample_name  ):
		sample_name = 'X'+sample_name
	return sample_name
def NameTransBack( sample_name ):
	if re.search( '(^X\d+)',sample_name  ):
		sample_name = sample_name[1:]
	sample_name = sample_name.replace("__"," vs ")
	return sample_name

if not os.path.exists(  outputpath ):
	os.makedirs( outputpath )
outputpath = os.path.abspath( outputpath  ) + os.sep
if not os.path.exists( input_path  ):
	print( 'ERROR!!! THE Input expression PATH doesn\'t exits!!!!!'  )
	sys.exit()
output_hash = {}

for a,b,c  in os.walk(  input_path  ):
	for f in c:
		if f.endswith(append):
			name = f.rsplit(".",1)[0]
			output_hash[sampleNameTrans(name)] = a+'/'+f

cluster_all_file = "\n".join(
    [ """%s=read.delim( "%s", header=TRUE, stringsAsFactors=TRUE ) """%(name,output_hash[name])    for name in  sorted(output_hash)   ]
)

cluster_cont = "compare_all <- list( "+",".join([ "%s$id"%(name)  for name in  sorted(output_hash) ])+")"
cluster_name = "names(compare_all)<-c("+",".join( """ "%s"  """%( NameTransBack(name) )for name in sorted(output_hash))+")"
if len(output_hash)<2:
	sys.exit()
Script="""
library("clusterProfiler")
%s
%s
%s

width = 6*length(compare_all)
gomap = read.table(file = "%s",header=FALSE)
names(gomap)<-c("entrezgene","go_accession")
buildGOmap(gomap,compress=FALSE)
dev.new()
pdf("%s/All_Compare_All.pdf",width=width)
ck <- compareCluster(geneCluster = compare_all, fun = "enrichGO",organism = "test", ont = "BP", readable = FALSE)
clusterProfiler::plot(ck)
dev.off()
"""%(
    cluster_all_file,
    cluster_cont,
    cluster_name,
    gomapping,
    outputpath
   )
SCRIPT = open( outputpath+'/Cluster.R','w'  )
SCRIPT.write(Script)

