#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/3
from lpp import *
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

parser.add_option("-s", "--STATIC", action="store", 
                  dest="static_path", 

                  help="static_path of QC")

parser.add_option("-i", "--INPUT", action="store", 
                  dest="input_path", 

                  help="input path of expression")


parser.add_option("-t", "--THRESHOLD", action="store", 
                  type='float',
                  dest="threshold", 

                  help="the threshold of padj to be considered as significant!!")

parser.add_option("-p", "--PVALUE",  
                  action="store_true", dest="pvale", default=False,  
                  help="use pval instead padj to do staistics test!")  

parser.add_option("-g", "--GOMAPPING",  
                  dest="gomapping", default=True,  
                  help="A file to store gene and GO mapping!") 



(options, args) = parser.parse_args()

threshold = options.threshold

gomapping = options.gomapping


if not os.path.exists(gomapping):
	print( "GOMAPPING file not exist!" )
	sys.exit()
else:
	gomapping = os.path.abspath(gomapping)

pvale = options.pvale
if pvale:
	measure='pval'
else:
	measure='padj'


outputpath   = options.outputpath

static_path   = os.path.abspath( options.static_path )+os.sep

input_path   = options.input_path
 
append = options.append
static_path = options.static_path

if not os.path.exists(  outputpath ):
	os.makedirs( outputpath )

outputpath = os.path.abspath( outputpath  ) + os.sep

def sampleNameTrans( sample_name ):
	if re.search( '(^\d+)',sample_name  ):
		sample_name = 'X'+sample_name
	return sample_name

if not os.path.exists( static_path  ):
	print( 'ERROR!!! THE Static PATH doesn\'t exits!!!!!'  )
	sys.exit()

if not os.path.exists( input_path  ):
	print( 'ERROR!!! THE Input expression PATH doesn\'t exits!!!!!'  )
	sys.exit()
	
	



input_path = os.path.abspath(  input_path )+os.sep

# To store the total depth of Sequencing 
size_factor = {}

for each_f in glob.glob(static_path +'*.stats'):
	STATIC = open ( each_f ,'rU' )
	STATIC.next()
	sample_name = sampleNameTrans( os.path.split(each_f)[-1].split('.')[0] )

	size_factor[ sample_name ] = STATIC.next().split()[2]
	
for each_matrix in glob.glob(  input_path+'*.'+append  ):
	stats_name = os.path.split(each_matrix)[-1].split('.')[0]
	
	sample_list = [x  for x in stats_name.split('__')]
	
	
	end_path = outputpath+stats_name
	
	if not os.path.exists(  end_path ):
		os.makedirs( end_path )
	end_path = os.path.abspath(  end_path )+os.sep
		
	matrix_output_path   = outputpath + os.path.split(  each_matrix  )[-1]+os.sep+ each_matrix+os.sep
	MATRIX = open( each_matrix,'rU'  )
	name_list = MATRIX.next()[:-1].split('\t')
	x_name = sampleNameTrans( name_list[1] )
	
	y_name = sampleNameTrans( name_list[2] )
	
	
	x_coverage = size_factor[ x_name ]
	
	y_coverage = size_factor[ y_name  ]
	
	matrix_abspath = each_matrix
	
	end_prefix = end_path + stats_name
	
	
	
	r_script = '''#!/usr/bin/Rscript
# functions

exampleFile = "%(matrix_abspath)s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE ) 

jpeg(file="%(end_prefix)s_dis.jpeg",width = 1050,height=480)
par(mfrow=c(1,2) )
plot(
 countsTable$%(x_name)s,
 countsTable$%(y_name)s,
 xlab='%(x_name)s',
 ylab='%(y_name)s',
 main='%(x_name)s vs%(y_name)s',
 pch=20, cex=.1, 
  )

library( DESeq )


rownames( countsTable ) <- countsTable$gene  
countsTable <- countsTable[ , -1 ]
conds <- c( "T", "N" ) 
cds <- newCountDataSet( countsTable, conds ) 
libsizes <- c(%(x_name)s=%(x_coverage)s, %(y_name)s=%(y_coverage)s)
sizeFactors(cds) <- libsizes  
cds <- estimateSizeFactors( cds ) 
cds <- estimateDispersions( cds,method='blind',fitType='local',sharingMode="fit-only" )   
res <- nbinomTest( cds, "T", "N" ) 

plotDE <- function( res ){
 plot(
 res$baseMean,
 res$log2FoldChange,
 xlab = 'baseMean',
 ylab = 'baseMean',
 main = '%(y_name)s vs %(x_name)s Diff',
 
 log="x", pch=20, cex=.1,
 col = ifelse( res$%(pvale)s < %(threshold)s, "red", "black" ) )
 }
 
plotDE( res )
dev.off()
res<-res[!is.na(res$%(pvale)s ),   ]
resSig <- res[ res$%(pvale)s  < %(threshold)s, ]  
up_resSIG <-resSig[resSig$log2FoldChange>0,]
down_resSIG <-resSig[resSig$log2FoldChange<0,]
write.table(resSig,row.names=FALSE,file='%(end_prefix)s.end',quote=FALSE,sep='\t') 
write.table(up_resSIG,row.names=FALSE,file='%(end_prefix)s_up.end',quote=FALSE,sep='\t') 
write.table(down_resSIG,row.names=FALSE,file='%(end_prefix)s_down.end',quote=FALSE,sep='\t') 



#GO分析并出图 
library("clusterProfiler")

Analysis_func<-function(Sig,up_Sig,down_Sig,output_prefix){
	run_analysis <- function(Diff,o_prefix){
		dev.new()
		GOe <- enrichGO(Diff$id, organism = "test", ont = "BP",pvalueCutoff = 0.05, qvalue = 0.1, readable = FALSE)
		if (nrow(summary(GOe))){
			t_prefix = paste(o_prefix,"_Enrich",sep="")
			pdf(paste(t_prefix,".pdf",sep=""))
			#postscript(paste(t_prefix,".eps",sep=""))
			#tiff( paste(t_prefix,".tiff",sep=""), res = 300)
			
			enrichMap(GOe)
			


			dev.off()
			#画concept图
			dev.new()
			t_prefix = paste(o_prefix,"_concept",sep="")
			pdf(paste(t_prefix,".pdf",sep=""))
			#postscript(paste(t_prefix,".eps",sep=""))
			#tiff( paste(t_prefix,".tiff",sep=""), res = 300)
			foldchange <-Diff$foldChange
			names(foldchange)<-Diff$id
			cnetplot(GOe, categorySize="geneNum",foldChange=foldchange)
			

			dev.off()
			write.table(summary(GOe), file = paste(o_prefix,"enrichment.tsv", sep = ""),sep="\t",quote=FALSE)
			}
		}
		
	run_analysis(Sig,paste(output_prefix ,"_ALL",sep="")   )
	run_analysis(up_Sig,paste(output_prefix ,"-Down",sep="")   )
	run_analysis(down_Sig,paste(output_prefix ,"-Up",sep="")   )
	
	compare_all <- list(  Sig$id,up_Sig$id,down_Sig$id)
	names(compare_all)<-c(paste(output_prefix,"_All",sep=""),paste(output_prefix,"-Down",sep=""),paste(output_prefix,"-Up",sep=""))
	dev.new()
	width <-length(compare_all)*10
	pdf(paste(output_prefix,"_compare.pdf",sep=""),width = width)
	#postscript(paste(output_prefix,"_compare.eps",sep=""))
	#tiff( paste(output_prefix,"_compare.tiff",sep=""), res = 300)
	ck <- compareCluster(geneCluster = compare_all, fun = "enrichGO",organism = "test", ont = "BP", readable = FALSE)
	plot(ck)
	dev.off()
	
}

gomap = read.table(file = "%(go_mapping)s",header=FALSE)
new_gomap[,1]<-gomap[,2]
new_gomap[,2]<-gomap[,1]
gomap<-new_gomap
buildGOmap(gomap,compress=FALSE)
Analysis_func(resSig,up_resSIG,down_resSIG,"%(end_prefix)s")










 '''%( 
	    {

	        "end_prefix":end_prefix,
	        "matrix_abspath":matrix_abspath,
	        "x_name":x_name,
	        "y_name":y_name,
	        "x_coverage":x_coverage,
	        "y_coverage":y_coverage,
	        "threshold":threshold,
	        "pvale":measure,
	        "go_mapping":gomapping,

	    } 
	)
	RSCRIPT = open( end_path+'Deseq.R'   ,'w' )
	RSCRIPT.write( r_script )
	RSCRIPT.close()
	#os.system(  '/usr/bin/Rscript ' +RSCRIPT.name+'&' )

	





