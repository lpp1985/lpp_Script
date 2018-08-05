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

#parser.add_option("-d", "--DATA", action="store", 
#                  dest="data_path", 
#
#                  help="fastq_path of Data")

parser.add_option("-i", "--INPUT", action="store", 
                  dest="input_path", 

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


#data_path   = os.path.abspath( options.data_path )+os.sep

input_path   = options.input_path

append = options.append
#data_path = options.data_path

if not os.path.exists(  outputpath ):
    os.makedirs( outputpath )

outputpath = os.path.abspath( outputpath  ) + os.sep

def sampleNameTrans( sample_name ):
    if re.search( '(^\d+)',sample_name  ):
        sample_name = 'X'+sample_name
    return sample_name

#if not os.path.exists( data_path  ):
#    print( 'ERROR!!! THE Static PATH doesn\'t exits!!!!!'  )
#    sys.exit()

if not os.path.exists( input_path  ):
    print( 'ERROR!!! THE Input expression PATH doesn\'t exits!!!!!'  )
    sys.exit()





input_path = os.path.abspath(  input_path )+os.sep

# To store the total depth of Sequencing 
#size_factor = {}

#for a,b,c in os.walk(data_path):
#    for e_f in c:
#        if e_f.endswith(".pair1"):
#            line_num = int(os.popen("wc -l %s"%(a+'/'+e_f)).read().split()[0])/2
#        
#            sample_name = sampleNameTrans( os.path.split(e_f)[-1].split('.')[0] )
#        
#            size_factor[ sample_name ] = str( line_num )
#print( size_factor )
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


    #x_coverage = size_factor[ x_name ]

    #y_coverage = size_factor[ y_name  ]

    matrix_abspath = each_matrix

    end_prefix = end_path + stats_name



    r_script = '''#!/usr/bin/Rscript
# functions
require(ggplot2)
require(ggthemes)
require(grid)

library( DESeq )
exampleFile = "%(matrix_abspath)s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE ) 

pdf(file="%(out_prefix)s_dis.pdf",width = 20)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

p<- ggplot(countsTable,aes(%(x_name)s,%(y_name)s))+geom_point()+xlab("%(x_name)s ReadCount")+theme_few()+ylab("%(y_name)s ReadCount")+ggtitle("%(x_name)s vs %(y_name)s")

vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}

print(p, vp = vplayout(1,1))

rownames( countsTable ) <- countsTable$gene  
countsTable <- countsTable[ , -1 ]
conds <- c( "T", "N" ) 
cds <- newCountDataSet( countsTable, conds ) 
libsizes <- c(%(x_name)s=sum(countsTable$%(x_name)s), %(y_name)s=sum(countsTable$%(y_name)s ) )
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

dev.new()
pdf(file="%(out_prefix)s_RPKM.pdf")

p3<- ggplot(res,aes(log10(baseMeanA),log10(baseMeanB),col=Condition))+geom_point()+ylab("%(y_name)s baseMean")+theme_few()+xlab("%(x_name)s baseMean")+scale_colour_manual(values=c("green", "blue", "red"))+ guides(colour = guide_legend(title = "FDR<0.05 and |log2Foldchange|>=1"))+theme(legend.position=c(.2, .9))+xlim(-2,max(log10(res$baseMean))+1)
ggplot_build(p3)
dev.off()

'''%(
       {
           "matrix_abspath":matrix_abspath,
           "out_prefix":end_prefix,
           "x_name":x_name,
           "y_name":y_name,
           "para":para,
           "threshold":threshold

       }


   )
    RSCRIPT = open( end_path+'Deseq.R'   ,'w' )
    RSCRIPT.write( r_script )
    RSCRIPT.close()
    os.system(  'Rscript ' +RSCRIPT.name +'&')







