#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/9
"""
from Taxon_GI_Parse import *
if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] 
         parse eggNOG data
   
         '''
    parser = OptionParser(usage =usage )    
    parser.add_option("-I", "--InputData", action="store",
                      dest="InputData",
                      help="blastn result of alignment")
    parser.add_option("-Y", "--Animal", action="store",
                      dest="Animal",
                      help="The Animal QueryId")
    parser.add_option("-o", "--Other", action="store",
                      dest="Other",
                      help="The Other QueryID")    
   
    
    (options, args) = parser.parse_args()
    
    RAW = open(options.InputData,'rU')
    animal_hash = {}
    other_hash = {}    
    for line in RAW:
        line_l = line.split("\t")
        name = line_l[2].split()[0]
        try:
            gi = int(re.search("gi\|(\d+)",line_l[5]).group(1))
        except:
            print(line_l[5])
        taxon_get = Taxon_GI.selectBy(GI=gi)

        Allanimal = Taxon_Classification.selectBy(Class="Animal")
        if taxon_get.count():
            taxon = taxon_get[0].Taxon
            has_animal =Taxon_Classification.selectBy(Class="Animal",Taxon = taxon) 
            if has_animal.count():
                animal_hash[name] = ""
            else:
                other_hash[name] = ""
        else:
            animal_hash[name] = ""
    
ANI = open(options.Animal,'w')
for key in animal_hash:
    ANI.write(key+'\n')
OTHER = open(options.Other,'w')
for key in other_hash:
    OTHER.write(key+'\n')
