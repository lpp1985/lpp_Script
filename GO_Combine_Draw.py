#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/5/4
"""

from lpp import *
from optparse import OptionParser
if __name__=='__main__':
    usage = '''usage: python2.7 %prog [options] Kmer
    



    Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
    parser = OptionParser(usage =usage )



    parser.add_option("-o", "--OUTPUT", action="store",
                      dest="output",
                      help="Output")  
    (options, args) = parser.parse_args()
    output = options.output

    END = open("%s.matrix"%(os.getpid()),'w')
    END.write("Situation\tTerm\tOntology\tGeneRatio\tQ_value\n")
    
    Outputhash_name = Ddict()
    sample_list = {}
    for f in args:
        RAW = pd.read_table(f)
    
        RAW["GeneRatio"] = RAW["numDEInCat"]/RAW["numInCat"]
    
        new_table = pd.DataFrame(RAW,columns=["term","ontology","GeneRatio","qvalue"])
        sample_name = os.path.dirname(f).rsplit("/",1)[-1]
        for i in xrange( 0, len( new_table )  ):
            data = new_table.iloc[i]
            try:
            	END.write(sample_name+'\t'+"\t".join(  [ data["term"],data["ontology"] ,str(data["GeneRatio"]),str(data["qvalue"] )] )+'\n'   )
            except:
		print( data['term'],f  )
     
    R_CACHE = open("%s.R"%(os.getpid()),'w')
    R_CACHE.write("""
    library(ggplot2)
    library(stringr)
    require(ggthemes)
    countsTable <- read.delim( "%(inp)s", header=TRUE, stringsAsFactors=TRUE ) 
    pathway_size = length(levels(factor(countsTable$Term)))
    Sample_size = length(levels(factor(countsTable$Situation)))
    if (pathway_size<10) pathway_size<-10
    countsTable$Term <-  str_wrap(countsTable$Term, width = 30)
    pdf("%(out)s",width=4*Sample_size,height = 1*pathway_size)
    qplot(data = countsTable,x=Situation,y=Term,size=GeneRatio,color=Q_value)+scale_colour_gradient(low="red", high="blue")+facet_grid(Ontology~.,scales="free_y",space="free")+theme_few()+theme(axis.text.y = element_text(angle = 350, vjust = .8), axis.text.x = element_text(angle = 270)  )
    dev.off()
    
    
    
    """%(
           {
               "inp" :END.name,
               "out": output
           }
          
       
       )
       
       
       )
    os.system("Rscript %s"%(R_CACHE.name))
        
    #os.remove(R_CACHE.name)
    #os.remove(END.name)
