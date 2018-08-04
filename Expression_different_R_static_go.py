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
                  action="store_true", dest="pvale", default=True,  
                  help="use pval instead padj to do staistics test!")  

parser.add_option("-g", "--GOMAPPING",  
                  dest="gomapping", default=True,  
                  help="A file to store gene and GO mapping!") 


parser.add_option("-n", "--Annotation",  
                  dest="annotation", default=True,  
                  help="A file to store gene annotation!") 


(options, args) = parser.parse_args()

threshold = options.threshold

gomapping = options.gomapping

annotation = options.annotation

print( annotation )
if not os.path.exists(gomapping):
	print( "GOMAPPING file not exist!" )
	sys.exit()
else:
	gomapping = os.path.abspath(gomapping)
print( gomapping )
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
	#STATIC.next()
	sample_name = sampleNameTrans( os.path.split(each_f)[-1].split('.')[0] )

	size_factor[ sample_name ] = STATIC.next()
	
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

countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE ) 

rownames( countsTable ) <- countsTable$gene  
countsTable <- countsTable[ , -1 ]
conds <- c( "T", "N" ) 
cds <- newCountDataSet( countsTable, conds ) 
libsizes <- c(%(x_name)s=%(x_coverage)s, %(y_name)s=%(y_coverage)s)
sizeFactors(cds) <- libsizes  
cds <- estimateSizeFactors( cds ) 
cds <- estimateDispersions( cds,method='blind',fitType='local',sharingMode="fit-only" )   
res <- nbinomTest( cds, "T", "N" ) 

plotDE <- function( res )
 plot(
 res$baseMean,
 res$log2FoldChange,
 xlab = 'baseMean',
 ylab = 'baseMean',
 main = '%(y_name)s vs %(x_name)s Diff',
 
 log="x", pch=20, cex=.1,
 col = ifelse( res$%(pvale)s < %(threshold)s, "red", "black" ) )
 
 plotDE( res )
 dev.off()
 res<-res[!is.na(res$%(pvale)s ),   ]
 resSig <- res[ res$%(pvale)s  < %(threshold)s, ]  

 write.table(resSig,row.names=FALSE,file='%(end_prefix)s.end',quote=FALSE,sep='\t') 
 
#GO检验
 library(topGO)
#读取所有表达的基因

allPvalue<-resSig[[dim( res  )[2]] ]
names(allPvalue)<- resSig[[1]]
#筛选差异基因
topDiffGenes <- function(allScore) {
    alldiff<-res[res$%(pvale)s<%(threshold)s,]
    cc <-names(allScore) %(tag)s names( allScore[ alldiff[[1]] ] )
	names(cc)<- names( allScore )
	cc[cc==TRUE]<-1
	cc[cc!=TRUE]<-0
	return(cc)
    }
	
upDiffGenes <- function(allScore) {
	upgene<-res[res$%(pvale)s<%(threshold)s & res$baseMeanA <res$baseMeanB,]
	
    cc <-  names(allScore) %(tag)s names(allScore[upgene[[1]]] )  
	cc[cc==TRUE]<-1
	cc[cc!=TRUE]<-0
	names(cc)<- names( allScore )
	return(cc)
    }

downDiffGenes <- function(allScore) {
   	downgene<-res[res$%(pvale)s<%(threshold)s & res$baseMeanA >res$baseMeanB,]
	
    cc <-    names( allScore ) %(tag)s names( allScore[downgene[[1]]]  )
	cc[cc==TRUE]<-1
	cc[cc!=TRUE]<-0
	names(cc)<- names( allScore )
	return(cc) 
    }
geneID2GO <- readMappings(file = "%(go_mapping)s")


GENEAnnotation <-read.table( "%(annotation)s", sep="\\t",header=TRUE,stringsAsFactors=FALSE)
dbGetQuery(GO_dbconn(), paste("SELECT go_id,term FROM go_term "))->go_all 


getdata <-function(allScore){
    return(allScore>0)
	}
	
upPvalue = upDiffGenes(  allPvalue )
downPvalue = downDiffGenes( allPvalue )
diffPvalue = topDiffGenes( allPvalue )
geneID2GO_diff <- geneID2GO[ names( geneID2GO ) %(tag)s names(diffPvalue[diffPvalue==1])    ]
geneID2GO_up <- geneID2GO[ names( geneID2GO ) %(tag)s names(upPvalue[upPvalue==1])    ]
geneID2GO_down <- geneID2GO[ names( geneID2GO ) %(tag)s names(downPvalue[downPvalue==1])    ]

GO2geneID <- inverseList(geneID2GO_diff)

GOdata <- new("topGOdata",
    description = "Diff", ontology = "BP",
    allGenes = diffPvalue, 
    geneSel = getdata ,
    nodeSize = 10,
   annot = annFUN.gene2GO, 
    gene2GO = geneID2GO
)

UPGOdata <- new("topGOdata",
    description = "Upgrade", ontology = "BP",
    allGenes = upPvalue, 
    geneSel = getdata ,
    nodeSize = 10,
   annot = annFUN.gene2GO, 
    gene2GO = geneID2GO
)

DownGOdata <- new("topGOdata",
    description = "Down", ontology = "BP",
    allGenes = downPvalue, 
    geneSel = getdata ,
    nodeSize = 10,
   annot = annFUN.gene2GO, 
    gene2GO = geneID2GO
)


#统计学检验
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")



allRes <- GenTable(GOdata, classicFisher = resultFisher,
 classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "classicFisher", ranksOf = "classicFisher")
resDATA<- allRes[ allRes$classicFisher <=.05 ,  ]




write.table(resDATA,row.names=FALSE,file='%(end_prefix)s.go',quote=FALSE,sep='\t')
result <-c()

all_need<-GO2geneID[names(GO2geneID) %(tag)s allRes["GO.ID"][[1]]]
for (i in 1:length(all_need) ){     
 for (j in 1:length( all_need[[i]] )  ){
 resneed<- res[res$id==all_need[[i]][j],]
 result<- rbind(result, c(names( all_need[i]), go_all[go_all[1]==names(all_need[i]),][2], all_need[[i]][j]  
 ,GENEAnnotation[GENEAnnotation$Gene==all_need[[i]][j],][[2]],resneed$baseMeanA,resneed$baseMeanB
 
 
 
 )  )  
 }
 }
 
#colnames(result)<-c("GO","Term","Gene","Annotation","BaseMeanA","BaseMeanB")

write.table(result,row.names=FALSE,file='%(end_prefix)s.go_detail',quote=FALSE,sep='\t')








resultFisher <- runTest(UPGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(UPGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(UPGOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(UPGOdata, classicFisher = resultFisher,
    classicKS = resultKS, elimKS = resultKS.elim,
	orderBy = "classicFisher", ranksOf = "classicFisher")
resDATA<- allRes[ allRes$classicFisher <=0.05 ,  ]
write.table(resDATA,row.names=FALSE,file='%(end_prefix)s.up_go',quote=FALSE,sep='\t')
GO2geneID <- inverseList(geneID2GO_up)
result<-c()
all_need<-GO2geneID[names(GO2geneID) %(tag)s allRes["GO.ID"][[1]]]
for (i in 1:length(all_need) ){     
 for (j in 1:length( all_need[[i]] )  ){
 resneed<- res[res$id==all_need[[i]][j],]
 result<- rbind(result, c(names( all_need[i]), go_all[go_all[1]==names(all_need[i]),][2], all_need[[i]][j]  
 ,GENEAnnotation[GENEAnnotation$Gene==all_need[[i]][j],][[2]],resneed$baseMeanA,resneed$baseMeanB
 
 
 
 )  )  
 }
 }
 
#colnames(result)<-c("GO","Term","Gene","Annotation","BaseMeanA","BaseMeanB")

write.table(result,row.names=FALSE,file='%(end_prefix)s.up_detail',quote=FALSE,sep='\t')




resultFisher <- runTest(DownGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(DownGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(DownGOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(DownGOdata, classicFisher = resultFisher,
 classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "classicFisher", ranksOf = "classicFisher")
resDATA<- allRes[ allRes$classicFisher <=0.05 ,  ]
write.table(resDATA,row.names=FALSE,file='%(end_prefix)s.down_go',quote=FALSE,sep='\t')


GO2geneID <- inverseList(geneID2GO_down)

all_need<-GO2geneID[names(GO2geneID) %(tag)s  resDATA["GO.ID"][[1]]]
result<-c()
for (i in 1:length(all_need) ){     
 for (j in 1:length( all_need[[i]] )  ){
 resneed<- res[res$id==all_need[[i]][j],]
 result<- rbind(result, c(names( all_need[i]), go_all[go_all[1]==names(all_need[i]),][2], all_need[[i]][j]  
 ,GENEAnnotation[GENEAnnotation$Gene==all_need[[i]][j],][[2]],resneed$baseMeanA,resneed$baseMeanB
 
 
 
 )  )  
 }
 }
 
#colnames(result)<-c("GO","Term","Gene","Annotation","BaseMeanA","BaseMeanB")

write.table(result,row.names=FALSE,file='%(end_prefix)s.down_detail',quote=FALSE,sep='\t')









 '''%( 
	    {
	        "tag":"%in%",
	        "end_prefix":end_prefix,
	        "matrix_abspath":matrix_abspath,
	        "x_name":x_name,
	        "y_name":y_name,
	        "x_coverage":x_coverage,
	        "y_coverage":y_coverage,
	        "threshold":threshold,
	        "pvale":measure,
	        "go_mapping":gomapping,
	        "annotation":annotation
	    } 
	)
	RSCRIPT = open( end_path+'Deseq.R'   ,'w' )
	RSCRIPT.write( r_script )
	RSCRIPT.close()
	os.system(  '/usr/bin/Rscript ' +RSCRIPT.name+'&' )

	





