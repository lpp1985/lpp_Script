#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/4
"""

import re
import tempfile


from optparse import OptionParser
import poster,urllib2
from poster.encode import multipart_encode  
from poster.streaminghttp import register_openers 

register_openers() 
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-i", "--Sequence", action="store",
                  dest="Sequence",

                  help="Genome Sequence in fasta format")
parser.add_option("-o", "--out", action="store",
                  dest="outputprefix",

                  help="oututprefix")




if __name__ == '__main__':
    (options, args) = parser.parse_args()

    outputprefix = options.outputprefix
  

    url = "http://bigsdb.pasteur.fr/perl/bigsdb/bigsdb.pl"

    values = {

        "fasta_upload":open(options.Sequence,'rb'),
        "sequence": "",
        'submit':"Submit", 'db':"pubmlst_klebsiella_seqdef_public",'page' :"plugin", 'name':"RuleQuery", 'ruleset':'Sequence_Type_determination'

    }	
    datagen, headers = poster.encode.multipart_encode(values)



  
  
    req = urllib2.Request(url, data = datagen, headers = headers)
    aa=True
    job_data = urllib2.urlopen(req).read()
    job_id = re.search("(/perl/bigsdb/bigsdb\.pl\?db\=pubmlst_klebsiella_seqdef_public\&page=job\&id=BIGSdb_\S+)\"", job_data).group(1)
    status = True
    while status:
        
        result_data = urllib2.urlopen("http://bigsdb.pasteur.fr" + job_id).read()
        if "finished" in result_data:
            status = False
    
    END = open(options.outputprefix + ".html", 'w')
    END.write(result_data)