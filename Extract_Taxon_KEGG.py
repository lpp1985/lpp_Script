#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/7
"""

from lpp import *
from optparse import OptionParser

usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Sequence", action="store",
                  dest="Sequence",

                  help="KEGG PEP FILE")
parser.add_option("-T", "--Taxon", action="store",
                  dest="Taxon",

                  help="Taxonomy File")
parser.add_option("-o", "--Output", action="store",
                  dest="Output",

                  help="OutputPath")
(options, args) = parser.parse_args()
Sequence = options.Sequence
OUTPUT   = options.Output
Taxon   = options.Taxon
all_taxon = File_dict(open(Taxon,'rU')).read(1,1)
END = open(OUTPUT,'w')
for t,s in fasta_check(open(Sequence,'rU')):
	title = re.search("^>([^\:]+)", t).group(1)
	if title in all_taxon:
		END.write(t+s)
