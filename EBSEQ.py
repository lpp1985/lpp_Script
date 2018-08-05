#!/usr/bin/python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/6/3
from lpp import *
#usage python2.7 expression_matrix_build.py     [ count.py's end__1_append_file_name  ]  [   matrix_end   ]
import glob,re,sys
from optparse import OptionParser 
usage = ''' To do EBSEQ of Reads mapping matrix !
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

(options, args) = parser.parse_args()

threshold = options.threshold

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

    sample_list = [x  for x in stats_name.split('___')]


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
require(ggplot2)
require(ggthemes)
require(grid)
library("EBSeq")
exampleFile = "%(matrix_abspath)s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE ,row.names=1) 

pdf(file="%(out_prefix)s_dis.pdf",width = 20)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

plot_data <- data.frame(data1=log10(countsTable$%(x_name)s+1) ,data2=log10(countsTable$%(y_name)s+1) )


lm.ab<-lm(plot_data$data1 ~plot_data$data2 )
math = as.data.frame(coefficients(lm.ab))
rsq<-summary(lm.ab)$r.squared
rsq_label <-paste("r^2:"  , rsq )
p<- ggplot(countsTable,aes(log10(%(x_name)s+1),log10(%(y_name)s+1)))+geom_point()+xlab("%(x_name)s ReadCount")+theme_few()+ylab("%(y_name)s ReadCount")+ggtitle("%(x_name)s vs %(y_name)s")+ geom_abline(intercept = math[1,],slope=math[2,],col="red")+annotate("text", x=2, y=8, parse=TRUE,label=rsq_label  )


vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

print(p, vp = vplayout(1,1))
sizes=MedianNorm(countsTable)
factors = as.factor( names(sizes) )
countsTable <- as.matrix(countsTable)
EBOut=EBTest(Data=countsTable, Conditions=factors,sizeFactors=sizes, maxround=5)
EBDERes=GetDEResults(EBOut, FDR= %(threshold)s )
GeneFC=as.data.frame(PostFC(EBOut)$PostFC)

colnames(GeneFC)<-"FoldChange"

GeneFC$Log2FC<-log2(GeneFC$FoldChange)
MEAN <- as.data.frame(EBOut$C1Mean)
colnames(MEAN)<-"BaseMeanB"
MEAN2 <-as.data.frame(EBOut$C2Mean)
colnames(MEAN2)<-"BaseMeanA"
MEAN<-cbind( MEAN,MEAN2 )
Res<-data.frame( ID=rownames(EBDERes$PPMat)  )
Res2 <- as.data.frame(EBDERes$PPMat)
Res2 <-subset(Res2, select = -PPDE )
colnames(Res2)<-"FDR"
Res<-cbind(Res, Res2) 
%(append)s
Res<-cbind( Res,MEAN) 
Res<-cbind(Res, GeneFC) 
Res<-Res[!is.na(Res$FDR), ]
Res$BaseMean<- (Res$BaseMeanA+Res$BaseMeanB)/2 
Res$Condition = cbind(rep("Not DEGs",nrow(Res)))
upRes<-Res[Res$Log2FC>=1 & Res$FDR < %(threshold)s  ,  ]
Res$Condition[Res$Log2FC>=1 & Res$FDR < %(threshold)s    ]<-paste("Up regulated gene :",dim(upRes)[1])

downRes<-Res[Res$Log2FC<=-1 & Res$FDR < %(threshold)s  ,  ]
Res$Condition[Res$Log2FC<=-1 & Res$FDR < %(threshold)s   ]<-paste("Down regulated gene :",dim(downRes)[1])

upRes<-Res[Res$Log2FC>=1 & Res$FDR < %(threshold)s  ,  ]
downRes<-Res[Res$Log2FC<=-1 & Res$FDR < %(threshold)s  ,  ]

resSig <- Res[Res$FDR <%(threshold)s  & abs(Res$Log2FC  )>=1, ] 
write.table(resSig,row.names=FALSE,file='%(out_prefix)s.end',quote=FALSE,sep='\t')
write.table(upRes,row.names=FALSE,file='%(out_prefix)s_up.end',quote=FALSE,sep='\t')
write.table(downRes,row.names=FALSE,file='%(out_prefix)s_down.end',quote=FALSE,sep='\t')
p2<- ggplot(Res,aes(log2(BaseMean),Log2FC,col=Condition))+geom_point()+ylab("log2FoldChange")+theme_few()+xlab("baseMean")+ggtitle("%(x_name)s vs %(y_name)s Diff")+scale_colour_manual(values=c("green", "blue", "red"))+ geom_hline(yintercept=0,col="red")+ guides(colour = guide_legend(title = "FDR<%(threshold)s and |log2Foldchange|>=1"))+theme(legend.position=c(.2, .9))+xlim(-2,max(log2(Res$BaseMean))+1)
print(p2, vp = vplayout(1,2))
dev.off()



pdf(file="%(out_prefix)s_RPKM.pdf")

ggplot(Res,aes(log10(BaseMeanA),log10(BaseMeanB),col=Condition))+geom_point()+ylab("%(y_name)s baseMean")+theme_few()+xlab("%(x_name)s baseMean")+scale_colour_manual(values=c("green", "blue", "red"))+ guides(colour = guide_legend(title = "FDR<0.05 and |log2Foldchange|>=1"))+theme(legend.position=c(.2, .9))+xlim(-2,max(log10(Res$BaseMeanA))+1)
dev.off()

'''%(
       {
           "matrix_abspath":matrix_abspath,
           "out_prefix":end_prefix,
           "x_name":x_name,
           "y_name":y_name,
           "append":"Res<-Res[ Res$ID%in%rownames(MEAN),   ]",
           "threshold":threshold

       }


   )
    RSCRIPT = open( end_path+'Deseq.R'   ,'w' )
    RSCRIPT.write( r_script )
    RSCRIPT.close()
    os.system(  'Rscript ' +RSCRIPT.name +'&')







