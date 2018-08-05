#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/11/4
"""
from optparse import OptionParser
from lpp import *
from parse_eggNog import *

if __name__=="__main__":

    '# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
    usage='''usage: python %prog [options]

	    multiproecssing blast '''
    parser = OptionParser(usage =usage )

    parser.add_option("-i", "--Input", action="store",
                      dest="path",
                      type='string',
                      help="Input File")		




    parser.add_option("-o", "--Output", action="store", 
                      dest="output",
                      help="Output File prefix")




    (options, args) = parser.parse_args()
    gene_nog = NOG_GENE.select(NOG_GENE.q.NOG.startswith("COG")   )
    
    INPUT_DATA= open(options.path,'rU')
    tax_name = os.path.split(options.path)[-1].split(".")[0]
    data_hash = Ddict()
    title = INPUT_DATA.next()
    title_list = title.strip().split("\t")
    
    for line in INPUT_DATA:
        data_list =line.strip().split("\t")
        name = tax_name+'_'+data_list[0]
        for i in xrange(  1,len(data_list[1:])  ):
            data_hash[name][title_list[i]]=data_list[i]
            
PEP = open(options.output+".pep",'w')
NUCL = open(options.output+".nuc",'w')
FUNC  = open(options.output+".function",'w')
for gene,detail_hash in data_hash.items():
    if detail_hash["Kind"] == "CDS":
        PEP.write('>'+gene+'\n'+re.sub("\W+","",detail_hash["Seq_Protein"])+'\n')
        NUCL.write('>'+gene+'\n'+detail_hash["Seq_Nucleotide"]+'\n')
        need_nog = NOG_CAT.select(NOG_CAT.q.NOG==detail_hash["EggNOG_NOG"])
        all_cog_cat = "".join( x.Cat for x in need_nog   )
        FUNC.write("%s\t%s\t%s\n"%(gene,detail_hash["EggNOG_NOG"]+all_cog_cat,detail_hash["Function"]))