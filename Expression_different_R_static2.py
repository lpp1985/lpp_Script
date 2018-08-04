#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/3
from lpp import *
import itertools
import os
import pandas as pd
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



parser.add_option("-i", "--INPUT", action="store", 
                  dest="input", 

                  help="total count matrix")

parser.add_option("-c", "--CONDITION", action="store", 
                  dest="condition", 

                  help="input path of expression")

parser.add_option("-p", "--PARAMATER", action="store", 
                  dest="para", 
                  default="padj",

                  help="the paramater to do consideration!! padj or pval")
parser.add_option("-t", "--THRESHOLD", action="store", 
                  type='float',
                  dest="threshold", 

                  help="the threshold of padj to be considered as significant!!")

(options, args) = parser.parse_args()

threshold = options.threshold

outputpath   = options.outputpath

para = options.para
if para not in ["padj","pval"]:
    raise IOError


condition   = options.condition



input_data = os.path.abspath( options.input )


if not os.path.exists(  outputpath ):
    os.makedirs( outputpath )

outputpath = os.path.abspath( outputpath  ) + os.sep

def sampleNameTrans( sample_name ):
    if re.search( '(^\d+)',sample_name  ):
        sample_name = 'X'+sample_name
    return sample_name

if not os.path.exists( input_data  ):
    print( 'ERROR!!! THE Reads Count Matrix doesn\'t exits!!!!!'  )
    sys.exit()


all_ReadsCountMatrix = pd.read_table( input_data )

all_sample = list(all_ReadsCountMatrix.columns)[1:]


if not condition:
    for sample_list in itertools.combinations( all_sample,2 ):
    
    
        end_path = outputpath+"/"+"___".join(sample_list)+'/'
        stats_name = "___".join(sample_list)
        if not os.path.exists(  end_path ):
            os.makedirs( end_path )
        end_path = os.path.abspath(  end_path )+os.sep

        CONDITION = open( end_path+'/condiion.tsv','w')
        CONDITION.write( "Sample\tCondition\n")
        CONDITION.write(sample_list[0]+'\tControl\n')
        CONDITION.write(sample_list[1]+'\tTreated\n')
        CONDITION.close()
        
        x_name = sample_list[0] 
    
        y_name =  sample_list[1] 


    
        end_prefix = end_path + stats_name
    
    
    
        r_script = '''#!/usr/bin/Rscript
# functions
require(ggplot2)
require(ggthemes)
require(grid)

library( DESeq )
condFile = "%(condition)s"
colTable <- read.delim( condFile, header=TRUE, row.names="Sample",stringsAsFactors=TRUE )

exampleFile = "%(matrix_abspath)s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE, row.names=1  ) 
countsTable<- countsTable[,rownames(colTable)]

pdf(file="%(out_prefix)s_dis.pdf",width = 20)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

p<- ggplot(countsTable,aes(\"%(x_name)s\",\"%(y_name)s\"))+geom_point()+xlab("%(x_name)s ReadCount")+theme_few()+ylab("%(y_name)s ReadCount")+ggtitle("%(x_name)s vs %(y_name)s")

vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

print(p, vp = vplayout(1,1))
conds <- c( "T", "N" ) 
cds <- newCountDataSet( countsTable, conds ) 
libsizes <- c(%(x_name)s=sum(countsTable$%(x_name)s), %(y_name)s= sum( countsTable$%(y_name)s ) )
sizeFactors(cds) <- libsizes  
cds <- estimateSizeFactors( cds ) 
cds <- estimateDispersions( cds,method='blind',sharingMode="fit-only" ,fitType="local" )   
res <- nbinomTest( cds, "T", "N" ) 
res$Condition = cbind(rep("Not DEGs",nrow(res)))
res$Condition[res$foldChange>1 & res$%(para)s < %(threshold)s   ]<-"Up regulated gene"
res$Condition[res$foldChange<1 & res$%(para)s < %(threshold)s  ]<-"Down regulated gene"
p2<- ggplot(res,aes(log2(baseMean),log2FoldChange,col=Condition))+geom_point()+ylab("log2FoldChange")+theme_few()+xlab("baseMean")+ggtitle("%(x_name)s vs %(y_name)s Diff")+scale_colour_manual(values=c("green", "blue", "red"))+ geom_hline(yintercept=0,col="red")+ guides(colour = guide_legend(title = "FDR<0.05 and |log2Foldchange|>=1"))+theme(legend.position=c(.2, .9))+xlim(-2,max(log2(res$baseMean))+1)
print(p2, vp = vplayout(1,2))
dev.off()


resSig <- res[res$%(para)s < %(threshold)s, ]
resSig <- resSig[!is.na(resSig$id), ]
resSig <- resSig[abs(resSig$log2FoldChange)>1, ]
upSig<- resSig[resSig$foldChange>1,]
downSig<- resSig[resSig$foldChange<1,]
write.table(resSig,row.names=FALSE,file='%(out_prefix)s.end',quote=FALSE,sep='\t')
write.table(upSig,row.names=FALSE,file='%(out_prefix)s_up.end',quote=FALSE,sep='\t')
write.table(downSig,row.names=FALSE,file='%(out_prefix)s_down.end',quote=FALSE,sep='\t')

pdf(file="%(out_prefix)s_RPKM.pdf")

p3<- ggplot(res,aes(log10(baseMeanA),log10(baseMeanB),col=Condition))+geom_point()+ylab("%(y_name)s baseMean")+theme_few()+xlab("%(x_name)s baseMean")+scale_colour_manual(values=c("green", "blue", "red"))+ guides(colour = guide_legend(title = "FDR<0.05 and |log2Foldchange|>=1"))+theme(legend.position=c(.2, .9))+xlim(-2,max(log10(res$baseMean))+1)
ggplot_build(p3)
dev.off()
    
    '''%(
           {
               "matrix_abspath":input_data,
               "out_prefix":end_prefix,
               "x_name":x_name,
               "y_name":y_name,
               "condition":CONDITION.name,
               "para":para,
               "threshold":threshold
    
           }
    
    
       )
        RSCRIPT = open( end_path+'Deseq.R'   ,'w' )
        RSCRIPT.write( r_script )
        RSCRIPT.close()
        os.system(  'Rscript ' +RSCRIPT.name )







