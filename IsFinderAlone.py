#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/4
"""
from lpp import *
import os
from bs4 import BeautifulSoup
import tempfile
from numpy import average
import pandas as pd
from optparse import OptionParser
import poster,time,urllib2,urllib
from poster.encode import multipart_encode  
from poster.streaminghttp import register_openers 
from copy import copy
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
    DATA = fasta_check( open(options.Sequence,'rU') )
    outputprefix = options.outputprefix
    #outputpath = check_path(  os.path.dirname( outputprefix   )   )
    #README = open(outputpath+'/Readme','w')
    #README.write(   
#"""
#使用在线网站ISFINDER（https://www-is.biotoul.fr/）预测IS序列，使用blastn+e value 1e-5作为参数
#包含以下结果：
#*.fa预测到的IS序列
#*.xls IS预测结果的表格
#*.stat IS预测结果的统计结果，用excel打开！！



    #""")

    url = "https://www-is.biotoul.fr/blast/ncbiIS.php"
    values = {
        "title":"",
        "prog":"blastn",
        "blast":"ok",
        "wordsize":11,
        "database":"ISfindernt",
        "seqfile":open(options.Sequence,'rb'),
        "seq":"",
        "expect": "1e-5"	,
        "gapcosts":"1 1"
    }	
    datagen, headers = poster.encode.multipart_encode(values) 


    sequence = re.sub('\s+','',DATA.next()[-1])


    is_stat = Ddict()
    #data = urllib.urlencode(values)
    req = urllib2.Request(url,datagen, headers)
    response = urllib2.urlopen(req)
    #try:
    uploadend = response.read()

    out_url = re.search("""(resultat.php\S+\"\>)""", uploadend).group(1)

    result = None
    while not  result:
        time.sleep(5)
        end_output = urllib.urlopen("https://www-is.biotoul.fr/blast/"+out_url).read()
        if "Query=" in end_output:
            result = end_output.split("</article>")[0]


    if result:
        
        ALN = open( outputprefix,'w'  )


        
        i=0
        data_list  = result.split("<b>Query=")[1:]

        for e_b in data_list:
            e_b = e_b.replace("</td>","\t</td>").replace("</th></tr>","\n").replace("</th>","\t</th>")
            data = BeautifulSoup(e_b,"html5lib")

            data = data.get_text() 

            block_list = data.split("\n\n",2)
            source_name,alignmentblock,blastblock = block_list

            is_detail = {}
            for each_isline in alignmentblock.split("\n")[2:]:
                isline_l = each_isline.split('\t')

                is_detail = {    
                    "Kind":"IS_Element",
                    "IS_Family":isline_l[1],
                    "IS_Origin":isline_l[3],
                    "Function":isline_l[0],



                }
                break

            is_finalResult = Ddict()
            is_statsis = Ddict()
            all_lignblock = Ddict()
            for eachblast_block in blastblock.split('>')[1:]:
                
                
                alignment_list = eachblast_block.split( " Score = " )

                startdata = alignment_list[ 0 ].strip()
                subject_name = startdata.split()[0]
                
                SubjLength = re.search(  "Length=(\d+)", startdata).group(1)
         
                for each_blastdetail in alignment_list[1:2]:
                    tag = 0
                    blast_stats  = each_blastdetail
                   
                    bitcore = re.search( "^(\S+)\s+bits\s+", blast_stats).group(1)
                    e_value  = re.search( "\s+Expect\s+\=\s+(\S+)", blast_stats).group(1)
                    if float(e_value)>1e-5 :
                        break

                    Identities = re.search( "\s+Identities\s+\=\s+([^\,]+)\,",blast_stats).group(1)

                    Gaps = re.search( "\s+Gaps\s+\=\s+([^\n]+)\n",blast_stats).group(1)
                    strand_detail  = re.search( "\s+Strand\=(\S+)",blast_stats).group(1)
                   
                    if "Minus" in strand_detail:
                        Strand = '-'
                    else:
                        Strand = "+"
                  
                 
                    is_detail["IS_Bitscore"] = bitcore
                    is_detail["IS_Identities"] = Identities
                    is_detail["IS_Evalue"] = e_value
                    is_detail["IS_SubjectLength"] = SubjLength
                    is_detail["IS_Gaps"] = Gaps
                    is_detail["Ref_Frame"] = Strand
                    is_detail["Seq_Nucleotide"] = sequence
                    is_detail["Seq_Nucl_Length"] = str( len(sequence) )

                    break
            result_list = [

                is_detail["Kind"],
                is_detail["Function"],

                is_detail["Ref_Frame"],
                is_detail["Seq_Nucl_Length"]    ,  
                is_detail["Seq_Nucleotide"],
                is_detail["IS_Family"],
                is_detail["Function"],

                is_detail["IS_Origin"],
                is_detail["IS_Bitscore"],
                is_detail["IS_Evalue"],
                is_detail["IS_Identities"],
                is_detail["IS_Gaps"],
                is_detail["IS_SubjectLength"]
               
            ]
            ALN.write( '\t'.join(["Kind","Function","Ref_Frame","Seq_Nucl_Length","Seq_Nucleotide","IS_Family","IS_Group","IS_Origin","IS_Bitscore","IS_Evalue","IS_Identities","IS_Gaps","IS_SubjectLength"])+'\n' )
            ALN.write( "\t".join( result_list  )+'\n'  )
            ALN.close()
          
