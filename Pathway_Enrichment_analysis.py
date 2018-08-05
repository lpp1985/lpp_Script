#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/3/3
"""

from lpp import *
from optparse import OptionParser

def enrichment_analysis(diff_gene,total_anno_gene,sample_size,sampled_diff):
    #超几何分布检验该通路是否显著变化


    r_script = open('/tmp/phyper.%s.R' % os.getpid() ,'w')
    output = '/tmp/phyper.%s.dat' % os.getpid()
    r_script.write( 'result<- 1- phyper(%s, %s, %s, %s)\n'%(sampled_diff-1,diff_gene,total_anno_gene-diff_gene,sample_size))
    r_script.write("write.table(result,row.names=FALSE,file='%s',quote=FALSE,sep='\t') \n"%(output))
    r_script.close()
    os.system("Rscript %s"%(r_script.name))
    output = open(output,'rU')
    p_value = float(output.read().split("\n")[1])
    
    os.remove(r_script.name)
    os.remove(output.name)
    return p_value

def fdr(p_value_list):    
    input_data = open('/tmp/pval.%s.txt' % os.getpid() ,'w')
    
    input_data.write("pval\n")
    for key in p_value_list:
        input_data.write("%s\n"%(key))
    input_data.close()
    r_script = open('/tmp/fdr_qval.%s.R' % os.getpid() ,'w')
    output ='/tmp/fdr_qval.%s.dat' % os.getpid()
    r_script.write( 
        """library(qvalue)         
pData<-read.delim( "%s", header=TRUE, stringsAsFactors=TRUE ) 

        """%(
               input_data.name,
              
             
           )
    
    
    )
    r_script.write(
    """
    pData$padj<- p.adjust(pData$pval, method="BH")
    pData$qvalue<- qvalue(pData$pval,lambda=0)$qvalues
    write.table(pData,row.names=F,file='%s',quote=FALSE,sep='\t')
    
    """ %(   
        output
        )
    )
    r_script.close()
    os.system("Rscript %s"%(r_script.name))
    p_adj = []
    q_qval = []
    output = open(output)
    output.next()
    for line in output:
        line_l = line.split("\t")
        p_adj.append( float(line_l[-2]) )
        q_qval.append( float(line_l[-1]) )
    os.remove(output.name)
    os.remove(r_script.name)
    os.remove(input_data.name)
    return p_adj,q_qval



if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] Kmer
        
    
    
    
    Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
    parser = OptionParser(usage =usage )



    parser.add_option("-a", "--PathwayAll", action="store",
                      dest="all",
                      help="all")
    parser.add_option("-d", "--DiffGene", action="store",
                      dest="diff",
                      help="Diff")
    parser.add_option("-e", "--Enrich", action="store",
                      dest="enrich",
                      help="Enrich")    
    (options, args) = parser.parse_args()
    Pathway_All = options.all
    