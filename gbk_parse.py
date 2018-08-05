#!/usr/bin/env python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/1/8

from Bio import SeqIO
from lpp import *
from Bio.SeqRecord import SeqRecord
import itertools

RECORD = SeqIO.parse( sys.argv[1],'genbank' )
GENE = open( sys.argv[2]+'.nuc','w' )
PRO = open( sys.argv[2]+'.pep','w' )
ANNO = open(sys.argv[2]+'.function','w')
from Bio import Seq
SeqRecord.format
for each_data in RECORD:
    
    seq = each_data.seq
    
    #seq = re.sub( '(\w{70})','\\1\n',seq  )
    for each_feature in each_data.features:
        if each_feature.type == 'CDS':
            nul_seq = each_feature.location.extract( seq ) 
#            print(each_feature)
            name    = re.sub('\s+','',str(each_feature.qualifiers['locus_tag'][0]))
#	    if name =="AKL27_RS26305":
#	    	print( each_feature )
#		print(each_feature.location)
            start =  int(each_feature.qualifiers["codon_start"][0])-1
            product = str( each_feature.qualifiers['product'][0] )
	    if "translation" not in each_feature.qualifiers:
		continue
            protein = re.sub('\s+','',each_feature.qualifiers['translation'][0])
	    
            nul = re.sub( '\s+','', str(nul_seq ) )[:len(protein)*3]
            if len(nul)!= len(protein)*3:
		print(name)
		print(start)
		print("%s\t%s"%(len(nul), len(protein)))
            GENE.write( '>%s\n%s\n'%( name, nul )  )
            PRO.write(   '>%s\n%s\n'%( name, protein  )   )
	    ANNO.write("%s\t-\t%s\n"%(name,product  )  )
		



