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


input_path   = options.input_path
 
append = options.append


if not os.path.exists(  outputpath ):
	os.makedirs( outputpath )

outputpath = os.path.abspath( outputpath  ) + os.sep

def sampleNameTrans( sample_name ):
	if re.search( '(^\d+)',sample_name  ):
		sample_name = 'X'+sample_name
	return sample_name



if not os.path.exists( input_path  ):
	print( 'ERROR!!! THE Input expression PATH doesn\'t exits!!!!!'  )
	sys.exit()
	
	



input_path = os.path.abspath(  input_path )+os.sep




	
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
	

	
	matrix_abspath = each_matrix
	
	end_prefix = end_path + stats_name
	
	
	
	r_script = '''#!/usr/bin/Rscript
# functions
library( DESeq2 )
options(bitmapType="cairo") 
exampleFile = "%(matrix_abspath)s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE ) 
x1<-matrix(data=countsTable$%(x_name)s,nc=1)
x2<-matrix(data=countsTable$%(y_name)s,nc=1)
countsTable2<-cbind(x1,x2)
dimnames(countsTable2)<- list(countsTable$gene,names(countsTable)[-1])


colData<-data.frame(condition=c("treated","untreated"),type=c("paired-end","paired-end"      ))
rownames(colData)<-c('%(x_name)s','%(y_name)s')
dds<-DESeqDataSetFromMatrix(countData= countsTable2,colData= colData,design =~condition)
dds<-DESeq(dds)
res<-results(dds)


pdf(file="%(end_prefix)s_dis.pdf")
par(mfrow=c(1,2) )
plot(
 countsTable$%(x_name)s,
 countsTable$%(y_name)s,
 xlab='%(x_name)s',
 ylab='%(y_name)s',
 main='%(x_name)s vs %(y_name)s',
 pch=20, cex=.1, 
  )
plotMA(res, main="%(x_name)s vs%(y_name)s  volcano Graph ", ylim=c(-2,2))
dev.off()
res_result <-as.data.frame(res) 
res_result<-res_result[!is.na(res$padj ),   ]
resSig <- res_result[ res_result$padj  < 0.1, ]
resSig$id<-row.names(resSig)
if (nrow(resSig)>0){

	up_resSIG <-resSig[resSig$log2FoldChange>0,]
	down_resSIG <-resSig[resSig$log2FoldChange<0,]
	write.table(resSig,row.names=FALSE,file='%(end_prefix)s.end',quote=FALSE,sep='\t') 
	write.table(up_resSIG,row.names=FALSE,file='%(end_prefix)s_up_end.tsv',quote=FALSE,sep='\t') 
	write.table(down_resSIG,row.names=FALSE,file='%(end_prefix)s_down_end.tsv',quote=FALSE,sep='\t') 



	#GO分析并出图 
	library("clusterProfiler")

	Analysis_func<-function(Sig,up_Sig,down_Sig,output_prefix){
		run_analysis <- function(Diff,o_prefix){
			
			GOe <- enrichGO(Diff$id, organism = "test", ont = "BP",pvalueCutoff = 0.5, qvalue = 0.5, readable = FALSE)
	        need_name = paste(o_prefix,"_EnrichBarplot",sep="")
	        
			if (!is.na(GOe)){

	            
	            

	            dev.new()
				t_prefix = paste(o_prefix,"_Enrich",sep="")
				pdf(paste(t_prefix,".pdf",sep=""))
				#postscript(paste(t_prefix,".eps",sep=""))
				#tiff( paste(t_prefix,".tiff",sep=""), res = 300)
				
				enrichMap(GOe,vertex.label.cex=.6)
				


				dev.off()
				#画concept图
				dev.new()
				t_prefix = paste(o_prefix,"_concept",sep="")
				pdf(paste(t_prefix,".pdf",sep=""))

				foldchange <-Diff$log2FoldChange
				head(foldchange)
				names(foldchange)<-Diff$id
				cnetplot(GOe, categorySize="geneNum",foldChange=foldchange,vertex.label.cex=.6)
				

				dev.off()
				write.table(summary(GOe),row.names=FALSE, file = paste(o_prefix,"enrichment.tsv", sep = ""),sep="\t",quote=FALSE)
				}
	        
			}
		
		run_analysis(Sig,paste(output_prefix ,"_ALL",sep="")   )

		run_analysis(up_Sig,paste(output_prefix ,"-Up",sep="")   )

	    run_analysis(down_Sig,paste(output_prefix ,"-Down",sep="")   )

		compare_all <- list(  Sig$id,up_Sig$id,down_Sig$id)
		
		ck_prefix = "%(name)s"
		names(compare_all)<-c(paste(ck_prefix,"_All",sep=""),paste(ck_prefix,"-Down",sep=""),paste(ck_prefix,"-Up",sep=""))
		return (  compare_all     )     
		
	}
	
	
	
	

	gomap = read.table(file = "%(go_mapping)s",header=FALSE)
	new_gomap <- gomap
	new_gomap[,1]<-gomap[,2]
	new_gomap[,2]<-gomap[,1]
	gomap<-new_gomap
	names(gomap)<-c("entrezgene","go_accession")
	buildGOmap(gomap)
	compare_all<-Analysis_func(resSig,up_resSIG,down_resSIG,"%(end_prefix)s")
	
	

}else{
q()
}

width = 6*length(compare_all)
dev.new()
pdf(paste("%(end_prefix)s","_compare.pdf",sep=""),width=width)


ck <- compareCluster(geneCluster = compare_all, fun = "enrichGO",organism = "test", ont = "BP", readable = FALSE)
clusterProfiler::plot(ck)
dev.off()



dev.new()
GOe <- enrichGO(resSig$id, organism = "test", ont = "BP",pvalueCutoff = 0.5, qvalue = 0.5, readable = FALSE)
name = paste("%(end_prefix)s_ALL" ,"_EnrichBarplot.pdf",sep="")
pdf(name)
barplot(GOe, drop=TRUE, showCategory=40)

dev.off()

dev.new()
GOe <- enrichGO(up_resSIG$id, organism = "test", ont = "BP",pvalueCutoff = 0.5, qvalue = 0.5, readable = FALSE)
name = paste("%(end_prefix)s-Up" ,"_EnrichBarplot.pdf",sep="")
pdf(name)
barplot(GOe, drop=TRUE, showCategory=40)

dev.off()

dev.new()
GOe <- enrichGO(down_resSIG$id, organism = "test", ont = "BP",pvalueCutoff = 0.5, qvalue = 0.5, readable = FALSE)
name = paste("%(end_prefix)s-Down" ,"_EnrichBarplot.pdf",sep="")
pdf(name)
barplot(GOe, drop=TRUE, showCategory=40)

dev.off()

 '''%( 
	    {

	        "end_prefix":end_prefix,
	        "matrix_abspath":matrix_abspath,
	        "x_name":x_name,
	        "y_name":y_name,

	        "threshold":threshold,
	        "pvale":measure,
	        "go_mapping":gomapping,
	        "name":stats_name
	    } 
	)
	RSCRIPT = open( end_path+'Deseq.R'   ,'w' )
	RSCRIPT.write( r_script )
	RSCRIPT.close()
	os.system(  '/pub/SOFTWARE/Other/R-3.2.0/bin/Rscript ' +RSCRIPT.name+'&& rm %s/Rplot*.pdf'%(os.getcwd()) )

	





