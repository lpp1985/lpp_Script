#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/9
"""
from Taxon_GI_Parse import *
from lpp import *
import os
if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] 
         parse eggNOG data
   
         '''
    parser = OptionParser(usage =usage )    
    parser.add_option("-I", "--InputData", action="store",
                      dest="InputData",
                      help="blastn result of alignment")
   
    
    (options, args) = parser.parse_args()
    
    RAW = fasta_check(open(options.InputData,'rU'))
    base_dir = os.path.split(  os.path.abspath(options.InputData)  )[0]+'/'

    Animal = open(base_dir+"Animal.fa",'w')
    Bacter = open(base_dir+"Bacter.fa",'w')
    Plants = open(base_dir+"Plant.fa",'w')
    
    for t,s in RAW:

        gi = int(re.search("gi\|(\d+)",t).group(1))
        taxon_get = Taxon_GI.selectBy(GI=gi)
        if taxon_get.count():
            taxon = taxon_get[0].Taxon
            taxon_class  =Taxon_Classification.selectBy(Taxon = taxon) 
            if taxon_class.count():
                classify = taxon_class[0].Class
                if classify =="Plants":
                    Plants.write(t+s)
                elif classify=="Bacter":
                    Bacter.write(t+s)
                elif classify =="Animal":
                    Animal.write(t+s)
    
#build index
os.system("ghostz db -i %s -o %s"%(options.InputData,base_dir+'nr&')     )
os.system("ghostz db -i %s -o %s"%(Animal.name,base_dir+'Animal&')     )
os.system("ghostz db -i %s -o %s"%(Bacter.name,base_dir+'Bacter&')     )
os.system("ghostz db -i %s -o %s"%(Plants.name,base_dir+'Plants&')     )