#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/16
"""
from lpp import *
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-o", "--Output", action="store",
                  dest="OutputPrefix",

                  help="OutputPrefix")
parser.add_option("-t", "--Threshold", action="store",
                  dest="Threshold",
                  type="float",

                  help="rpkm threshold")

parser.add_option("-s", "--Fasta", action="store",
                  dest="Seq",

                  help="Sequence File")

parser.add_option("-r", "--Rpkm", action="store",
                  dest="Rpkm",

                  help="Rpkm File")


parser.add_option("-c", "--Count", action="store",
                  dest="Count",

                  help="Reads Count File")
if __name__ == '__main__':
	(options, args) = parser.parse_args()
	thrshold = options.Threshold
	rpkm_data= pd.read_table(options.Rpkm)
	rpkm_data["max"] = rpkm_data[ rpkm_data.columns[1:]   ].max(1)
	rpkm_filter = rpkm_data[  rpkm_data["max"] >=thrshold   ]
	del rpkm_filter["max"]
	
	rpkm_filter.to_csv(options.OutputPrefix+'.rpkm',index = False,sep = "\t"   )
	
	
	########Draw Graph#######################
	RAW = open(options.OutputPrefix+'.rpkm','rU')
	RPDRAW = open( "R.Data",'w'  )
	RPDRAW.write("RPKM\tSample\n")
	title_l = RAW.next().strip().split("\t")[1:]
	
	for line in RAW:
		line_l = line.strip().split("\t")
		data_l = line_l[1:]
		for i in xrange(0,len(data_l)):
			RPDRAW.write(data_l[i]+'\t'+title_l[i]+'\n')
			
	RPDRAW.close()		
	R_Script = open("Draw.R",'w')
	R_Script.write(
	    """
library("ggplot2")
countsTable <- read.delim( "%s", header=TRUE, stringsAsFactors=TRUE )
pdf("%s_rpkm.pdf")
ggplot(countsTable, aes(log10(RPKM+1),fill=Sample))+geom_density(alpha=0.2) 
dev.off()


	    
"""%(
	   RPDRAW.name,
	   options.OutputPrefix
   
   
   )
	
	
	
	)
	os.system("Rscript %s"%(R_Script.name) )
	old_name = rpkm_filter.columns[1:]
	new_name = [ "RPKM_"+x for x in old_name    ]
	changname_hash = dict(zip(old_name,new_name))
	rpkm_filter.rename( columns= changname_hash, inplace=True  )
	
	all_filteredGene = list(rpkm_filter["Gene"])

	fil_geneHash =  dict(zip(all_filteredGene,[""]*len(all_filteredGene)))

	count_data = pd.read_table( options.Count )
	count_has_data = count_data[  count_data["Gene"].isin(fil_geneHash) ]	

	count_has_data.to_csv( options.OutputPrefix+'.count',sep="\t",index=False  )
	old_name = count_has_data.columns[1:]
	new_name = [ "ReadCount_"+x for x in old_name    ]
	changname_hash = dict( zip(old_name,new_name) )
	count_has_data.rename(columns= changname_hash, inplace=True  )	
	
	TMP = open("%s.tmp"%os.getpid(),'w')
	TMP.write("Gene\tSequence\n")
	SEQ = open(   options.OutputPrefix+'.fasta','w'    )
	BED = open(   options.OutputPrefix+'.bed','w'    )
	LENGTH = open(   options.OutputPrefix+'.length','w'    )
	for t,s in fasta_check(  open( options.Seq  )  ):
		name  = t[1:].split()[0]
		s1 = re.sub( "\s+", '', s )
		if name in fil_geneHash:
			LENGTH.write(name+'\t%s\n'%( len(s1) )  )
			TMP.write(name+'\t'+s1+'\n')
			SEQ.write('>'+name+'\n'+s)
			BED.write(
			 '\t'.join(
			    [
			        name,
			        "0",
			        str( len(s1)  ),
			        name,
			        '0',
			        '+',
			        '0',
			        str( len(s1)  ),
			        '0',
			        '1',
			        str( len(s1)+1  ),
			        '0'
			        
			    
			    
			    
			    ]
			 
			 
			 
			 )   +'\n'
			
			
			
			)
	TMP.close()
	seq_data = pd.read_table( TMP.name  )
	tsv_data = pd.DataFrame.merge( seq_data,rpkm_filter,on="Gene",how="left"   )
	tsv_data = pd.DataFrame.merge( tsv_data,count_has_data,on="Gene",how="left"   )
	tsv_data.rename(columns={tsv_data.columns[0]:"Name"},inplace = True)
	tsv_data.to_csv( options.OutputPrefix+'.xls',sep='\t',index=False)
	os.remove(TMP.name)
