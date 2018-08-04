#!/usr/bin/env python
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
parser.add_option("-g", "--GROUP", action="store", 
                  dest="group", 

                  help="the threshold of padj to be considered as significant!!")


(options, args) = parser.parse_args()

def sampleNameTrans( sample_name ):
	if re.search( '(^\d+$)',sample_name  ):
		sample_name = 'X'+sample_name
	return sample_name


group = options.group
threshold = options.threshold

outputpath   = options.outputpath

static_path   = os.path.abspath( options.static_path )+os.sep

input_path   = options.input_path
 
append = options.append


static_path = options.static_path

group_info = dict()
for line in open(group,'rU'):
	line_l = line[:-1].split()
	data_list = [ sampleNameTrans( x.strip() ) for x in line_l[-1].split(';')  ]
	for each_sample in data_list:
		group_info[ each_sample ] = line_l[0]






if not os.path.exists(  outputpath ):
	os.makedirs( outputpath )

outputpath = os.path.abspath( outputpath  ) + os.sep



if not os.path.exists( static_path  ):
	print( 'ERROR!!! THE Static PATH doesn\'t exits!!!!!'  )
	sys.exit()

if not os.path.exists( input_path  ):
	print( 'ERROR!!! THE Input expression PATH doesn\'t exits!!!!!'  )
	sys.exit()
	
	



input_path = os.path.abspath(  input_path )+os.sep

# To store the total depth of Sequencing 
size_factor = {}




for each_f in glob.glob(static_path +'*.total.stats'):
	STATIC = open ( each_f ,'rU' )
	STATIC.next()
	sample_name = sampleNameTrans( os.path.split(each_f)[-1].split('.')[0] )

	size_factor[ sample_name ] = STATIC.next().split('\t')[2]
	
for each_matrix in glob.glob(  input_path+'*.'+append  ):
	
	sample_mapping = {}
	
	already_have_group = {}
	
	seperate_status = []
	
	stats_name = os.path.split(each_matrix)[-1].split('.')[0]
	
	
	
	end_path = outputpath+stats_name
	
	
	
	if not os.path.exists(  end_path ):
		os.makedirs( end_path )
	end_path = os.path.abspath(  end_path )+os.sep
		
	matrix_output_path   = outputpath + os.path.split(  each_matrix  )[-1]+os.sep+ each_matrix+os.sep
	MATRIX = open( each_matrix,'rU'  )
	name_list =  [ sampleNameTrans(x) for x in      MATRIX.next()[:-1].split('\t')[1:]]

	conds = []
	libary_size = []
	already_have_group[   group_info[ name_list[0] ]    ] = ''
	for each_sp in name_list:
		each_sp = sampleNameTrans( each_sp )
		
		
		libary_size.append( '%s=%s'%( each_sp, size_factor[  each_sp ]  )  )
		
		if group_info[ each_sp ]  in  already_have_group:
			
			conds.append(  '"T"' )
			

		else:
			conds.append(  '"N"' )
	method  = 'pooled'
	for each_status in set( conds ):
		
		if conds.count( each_status ) ==1:
			method  = 'normal'
	conds = ','.join(  conds )
	
	libary_size = ',' .join(  libary_size  )
	
	matrix_abspath = each_matrix
	
	end_prefix = end_path + stats_name
	
	
	
	r_script = '''#!/usr/bin/Rscript
# functions

exampleFile = "%s"
countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE ) 

jpeg( file="%s_dis.jpeg" )



library( DESeq )

countsTable <- read.delim( exampleFile, header=TRUE, stringsAsFactors=TRUE ) 

rownames( countsTable ) <- countsTable$gene  
countsTable <- countsTable[ , -1 ]
conds <- c( %s ) 
cds <- newCountDataSet( countsTable, conds ) 
libsizes <- c(%s)
sizeFactors(cds) <- libsizes  
cds <- estimateSizeFactors( cds ) 
cds <- estimateVarianceFunctions( cds,method='%s' )  
res <- nbinomTest( cds, "T","N" ) 

plotDE <- function( res )
 plot(
 res$baseMean,
 res$log2FoldChange,
 xlab = 'baseMean',
 ylab = 'logbaseMean',

 
 log="x", pch=20, cex=.1,
 col = ifelse( res$padj < %s, "red", "black" ) )
 
 plotDE( res )
 dev.off()
 resSig <- res[ res$padj < %s, ]  

 write.table(resSig,row.names=FALSE,file='%s.end',quote=FALSE,sep='\t') '''%( matrix_abspath,end_prefix,conds, libary_size ,method,threshold,threshold, end_prefix )
	RSCRIPT = open( end_path+'Deseq.R'   ,'w' )
	RSCRIPT.write( r_script )
	RSCRIPT.close()
	#os.system(  RSCRIPT.name +'&')
	
	





