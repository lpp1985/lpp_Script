#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/4
"""
import ssl
from lpp import *
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
    outputpath = check_path(  os.path.dirname( outputprefix   )   )
    README = open(outputpath+'/Readme','w')
    README.write(   
"""
使用在线网站ISFINDER（https://www-is.biotoul.fr/）预测IS序列，使用blastn+e value 1e-5作为参数
包含以下结果：
*.fa预测到的IS序列
*.xls IS预测结果的表格
*.stat IS预测结果的统计结果，用excel打开！！



    """)

    url = "https://www-is.biotoul.fr/blast/ncbiIS.php"
    values = {
        "title":"",
        "prog":"blastn",
        "blast":"ok",
        "wordsize":11,
        "database":"ISfindernt",
        "seqfile":open(options.Sequence,'rb'),
        "seq":"",
        "expect": "1e-100"	,
        "gapcosts":"1 1"
    }	
    datagen, headers = poster.encode.multipart_encode(values) 


    sequence = re.sub('\s+','',DATA.next()[-1])
    NUL = open( outputprefix+".fa",'w'  )
    STAT = open( outputprefix+".stat",'w'  )
    HTML = open(outputprefix+".html",'w')
    ssl._create_default_https_context = ssl._create_unverified_context
    is_stat = Ddict()
    #data = urllib.urlencode(values)
    req = urllib2.Request(url,datagen, headers)
    aa=True
    while aa:
        try:
            response = urllib2.urlopen(req)
            aa = False
        except:
            pass
    #try:
    uploadend = response.read()

    out_url = re.search("""(resultat.php\S+\"\>)""", uploadend).group(1)

    result = None
    while not  result:
        time.sleep(5)
        while True:
            try:
                
                end_output = urllib.urlopen("https://www-is.biotoul.fr/blast/"+out_url).read()
                break
            except:
                pass
        if "Query=" in end_output:
            result = end_output.split("</article>")[0]
            HTML.write(end_output)

    if result:

        ALN = open( outputprefix+".xls",'w'  )
        STAT.write("IS_name\tNumber\tAverage.Length\n")

        ALN.write( '\t'.join(["Name","Ref_Source","Kind","Function","Ref_Start","Ref_Stop","Ref_Frame","Seq_Nucl_Length","Seq_Nucleotide","IS_Family","IS_Group","IS_Origin","IS_Bitscore","IS_Evalue","IS_Identities","IS_Gaps","IS_SubjectLength","IS_Frame"])+'\n' )
        i=0
        data_list  = result.split("<b>Query=")[1:]

        for e_b in data_list:
            e_b = e_b.replace("</td>","\t</td>").replace("</th></tr>","\n").replace("</th>","\t</th>")
            data = BeautifulSoup(e_b,"html5lib")

            data = data.get_text() 

            block_list = data.split("\n\n",2)
            source_name,alignmentblock,blastblock = block_list
            source_name = source_name.split()[0]
            is_detail = {}
            for each_isline in alignmentblock.split("\n")[2:]:
                isline_l = each_isline.split('\t')

                is_detail[ isline_l[0] ] = {    
                    "Kind":"IS_Element",
                    "IS_Group":isline_l[2],
                    "IS_Family":isline_l[1],
                    "IS_Origin":isline_l[3],
                    "Function":isline_l[0],
                    "Ref_Source":source_name



                }
            is_finalResult = Ddict()
            is_statsis = Ddict()
            all_lignblock = Ddict()
            for eachblast_block in blastblock.split('>')[1:]:


                alignment_list = eachblast_block.split( " Score = " )

                startdata = alignment_list[ 0 ].strip()
                subject_name = startdata.split()[0]

                SubjLength = re.search(  "Length=(\d+)", startdata).group(1)
                for each_blastdetail in alignment_list[1:]:
                    tag = 0
                    blast_stats , blast_detail = each_blastdetail.split("\n\n",1 )
                    bitcore = re.search( "\s+bits\s+\((\d+)\)", blast_stats).group(1)
                    e_value  = re.search( "\s+Expect\s+\=\s+(\S+)", blast_stats).group(1)
                    if float(e_value)>1e-5 :
                        continue

                    Identities = re.search( "\s+Identities\s+\=\s+([^\,]+)\,",blast_stats).group(1)

                    Gaps = re.search( "\s+Gaps\s+\=\s+([^\,]+)\n",blast_stats).group(1)
                    strand_detail  = re.search( "\s+Strand\=(\S+)",blast_stats).group(1)

                    if "Minus" in strand_detail:
                        Strand = '-'
                    else:
                        Strand = "+"
                    query_start = re.search("Query\s+(\d+)\s+", blast_detail).group(1)
                    query_end = re.findall("(\d+)\n+", blast_detail)[-2]
                    if not len(all_lignblock):
                        all_lignblock[int(query_start)][ int( query_end )]=""
                    else:

                        for start in all_lignblock:
                            if int(query_start )>=start:
                                for end in all_lignblock[start]:
                                    if end >= int( query_end ):
                                        tag = 1
                    if tag ==0:
                        all_lignblock[int(query_start)][ int( query_end )]=""
                    else:

                        continue
                    q_length =  abs(int(query_end) - int( query_start ))
                    if q_length > 6000:
                        continue
                    is_finalResult[int(query_start)][ subject_name ] = copy(is_detail[ subject_name ])
                    is_finalResult[int(query_start)][ subject_name ]["IS_Bitscore"] = bitcore
                    is_finalResult[int(query_start)][ subject_name ]["IS_Identities"] = Identities
                    is_finalResult[int(query_start)][ subject_name ]["IS_Evalue"] = e_value
                    is_finalResult[int(query_start)][ subject_name ]["IS_SubjectLength"] = SubjLength
                    is_finalResult[int(query_start)][ subject_name ]["IS_Gaps"] = Gaps
                    is_finalResult[int(query_start)][ subject_name ]["Ref_Frame"] = Strand
                    is_finalResult[int(query_start)][ subject_name ]["Seq_Nucl_Length"] = str( abs(int(query_end) - int( query_start ) ) )
                    is_finalResult[int(query_start)][ subject_name ]["Seq_Nucleotide"] = sequence[ int( query_start ) :int(query_end) ]
                    is_statsis[ subject_name  ][int(query_start)] = int(query_end) - int( query_start )
                    is_finalResult[int(query_start)][ subject_name ]["Ref_Start"] = query_start
                    is_finalResult[int(query_start)][ subject_name ]["Ref_Stop"] = query_end
            i=0

            for each_loc in sorted(is_finalResult):
                for each_result in sorted(is_finalResult[ each_loc  ] ):

                    i+=1
                    is_name = source_name+"_IS%s"%(i)

                    result_list = [
                        is_name,
                        is_finalResult[each_loc][each_result]["Ref_Source"],
                        is_finalResult[each_loc][each_result]["Kind"],
                        is_finalResult[each_loc][each_result]["Function"],
                        is_finalResult[each_loc][each_result]["Ref_Start"],
                        is_finalResult[each_loc][each_result]["Ref_Stop"],
                        is_finalResult[each_loc][each_result]["Ref_Frame"],
                        is_finalResult[each_loc][each_result]["Seq_Nucl_Length"],
                        is_finalResult[each_loc][each_result]["Seq_Nucleotide"],
                        is_finalResult[each_loc][each_result]["IS_Family"],
                        is_finalResult[each_loc][each_result]["IS_Group"],
                        is_finalResult[each_loc][each_result]["IS_Origin"],
                        is_finalResult[each_loc][each_result]["IS_Bitscore"],
                        is_finalResult[each_loc][each_result]["IS_Evalue"],
                        is_finalResult[each_loc][each_result]["IS_Identities"],
                        is_finalResult[each_loc][each_result]["IS_Gaps"],
                        is_finalResult[each_loc][each_result]["IS_SubjectLength"],



                    ]
                    ALN.write( "\t".join( result_list  )+'\n'  )
                    NUL.write('>'+is_name+'\n'+is_finalResult[each_loc][each_result]["Seq_Nucleotide"]+'\n')

            for key in is_statsis:
                length_all = []
                for each_element in is_statsis[key]:
                    length_all.append(is_statsis[ key ][ each_element ])
                STAT.write( 
                    "%s\t%s\t%s\n"% (  
                        key,len( is_statsis[key] ) , average( length_all ) 
                    )  
                )
            if not is_statsis:
                STAT.write(  "Not Find IS!!\n"  )	
    #except Exception,error:
        #print(error)

