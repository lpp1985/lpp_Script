#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/11/25
"""

#from lpp import *
import os,re
from collections import defaultdict
class Ddict(defaultdict,dict):
    def __init__(self):
        defaultdict.__init__(self, Ddict)
    def __repr__(self):
        return dict.__repr__(self)
class fasta_check(object):
    def __init__(self,file_handle):
        assert isinstance( file_handle,file ),'The Input paramater must be a File Handle'
        self.file=iter(file_handle)
        for line in self.file:
            if line[0]=='>':
                self.define=line
                break
    def __iter__(self):
        return self
    def next(self):
        if not self.define:
            raise StopIteration

        name=self.define
        self.define=''        
        s=[]
        for line in self.file:
            if line[0]!='>':
                s.append(line)
            else:
                self.define=line
                break
        s=''.join(s)
        return (name,s)

import pandas
from optparse import OptionParser
def Pagan( input_file  ):
    global Ortholog_Pair
    input_file = os.path.abspath( input_file  )
    path = os.path.split(input_file)[0]+'/'
    output = path+'outfile'
    ancestor = path+'ancestor'
    os.system("pagan --seq %s --threads 64 --silent -o %s   "  %(  input_file, output ))
    RAW = fasta_check( open( output+".fas" ,'rU'   ) )
    t,s = RAW.next()
    s = re.sub("\s+", "" , s)
    if s.count("-")/len(s)>0.4:
        Ortholog_Pair[orthId][ "AncestorySeq" ] = "-"
        return "Too Divergent!!"
    os.system("pagan --ref-seqfile %s.fas  --ref-treefile %s.tre  --output-ancestors -o %s  --threads 64 --silent"%( output,output,ancestor  ))
    for t,s in fasta_check(  open( "%s.fas"%(ancestor),'rU'   )  ):
        if "#1#" in t:
            s = re.sub("\s+", "", s)
            return s
        
def CdsFinder( input_file   ):
    Ortholog_Pair[orthId]["AncestorCDS"] = '-'
    Ortholog_Pair[orthId]["AncestorCDS_start"] = '-'
    Ortholog_Pair[orthId]["AncestorCDS_end"] = '-'
    Ortholog_Pair[orthId]["AncestorCDS_frame"] = '-'
    Ortholog_Pair[orthId]["mianCDS"] = '-'
    Ortholog_Pair[orthId]["mianCDS_end"] = '-'
    Ortholog_Pair[orthId]["mianCDS_frame"] = '-'  
    Ortholog_Pair[orthId]["mianCDS_start"] = '-'
    Ortholog_Pair[orthId]["yanCDS"] = '-'
    Ortholog_Pair[orthId]["yanCDS_start"] = '-'
    Ortholog_Pair[orthId]["yanCDS_end"] = '-'
    Ortholog_Pair[orthId]["yanCDS_frame"] = '-'
    
    
    path = os.path.split(input_file)[0]+'/'
    data_hash = Ddict()
    os.system(" TransDecoder -t %s 2>/dev/null 1>/dev/null --workdir  %s/prediction"%(input_file, path))
    size = os.path.getsize(  "%s.transdecoder.cds"%(input_file) )
    if size >0:
        CDS = fasta_check(open( "%s.transdecoder.cds"%(input_file),'rU'  ))
        for t,s in CDS:
            s = re.sub("\s+", "", s)
            name = t[1:].split("|")[0]
            [(start,end,frame)] = re.findall( "\:(\d+)\-(\d+)\((\S)\)",t  )
            if name not in data_hash:
                if s[-3:] in ["TAA","TAG","TGA"]:
                    s = s[:-3]
                data_hash[ name ]["Seq"] = s
                data_hash[ name ]["start"]=start
                data_hash[ name ]["end"]=end
                data_hash[ name ]["frame"]=frame
            elif len(s) >len(data_hash[ name ]["Seq"]):
                data_hash[ name ]["Seq"] = s
                data_hash[ name ]["start"]=start
                data_hash[ name ]["end"]=end
                data_hash[ name ]["frame"]=frame                
    if len(data_hash)==3:
        END = open(path+"all_cds.fa",'w')
        LOC = open(path+"Location.tsv",'w')
        for key in data_hash:
            if "Ancestor" not in key:
                ALL_SEQ.write('>'+key+'|CDS'+'\n'+data_hash[key][ "Seq" ])
                Ortholog_Pair[orthId]["AncestorCDS"] = data_hash[key][ "Seq" ]
                Ortholog_Pair[orthId]["AncestorCDS_start"] = data_hash[key][ "start" ]
                Ortholog_Pair[orthId]["AncestorCDS_end"] = data_hash[key][ "end" ]
                Ortholog_Pair[orthId]["AncestorCDS_frame"] = data_hash[key][ "frame" ]
            elif "mian" not in key:
                Ortholog_Pair[orthId]["mianCDS"] = data_hash[key][ "Seq" ]
                Ortholog_Pair[orthId]["mianCDS_start"] = data_hash[key][ "start" ]
                Ortholog_Pair[orthId]["mianCDS_end"] = data_hash[key][ "end" ]
                Ortholog_Pair[orthId]["mianCDS_frame"] = data_hash[key][ "frame" ]     
            else:
                Ortholog_Pair[orthId]["yanCDS"] = data_hash[key][ "Seq" ]
                Ortholog_Pair[orthId]["yanCDS_start"] = data_hash[key][ "start" ]
                Ortholog_Pair[orthId]["yanCDS_end"] = data_hash[key][ "end" ]
                Ortholog_Pair[orthId]["yanCDS_frame"] = data_hash[key][ "frame" ]                     
            END.write('>'+key+'\n'+data_hash[key][ "Seq" ]+'\n')
            
            LOC.write( key+'\t'+data_hash[key][ "start"  ]+'\t'+ data_hash[key][ "end"  ]+'\t'+ data_hash[key][ "frame"  ]+'\n')
        return END.name
    else:
        return "Not Enough CDS"
    
def KaksCal(  input_file  ):
    global Ortholog_Pair
    
    path = os.path.split(input_file)[0]+'/'
    output = path+"multiple_alignment"
    os.system("pagan --seq %s --threads 64 --silent -o %s   "  %(  input_file, output ))
    output_trimed = output+"_trimed.fas"
    os.system("""trimal -in  %s  -fasta |sed -r "s/\s+[0-9]+\s+bp//g" >%s """%( output+'.fas' ,output_trimed ) )
    cache_hash = {}
    for t,s in fasta_check( open( output_trimed,'rU'  )  ):
        s = re.sub("\s+", '', s)
        if "mian" in t or "M_" in t or "HARM_" in t:
            cache_hash["mian"] = s
        elif "yan" in t or "M_" in t or "HAS_" in t:
            cache_hash["yan"] = s
        else:
            cache_hash["Ances"] = s
    mian_ances_name = path+"mian_vs_ances.axt"
    M_A = open( mian_ances_name,'w' )
    M_A.write("Mian_vs_Ancestor\n")
    all_filter = []
    for loc_iter in re.finditer("\-+", cache_hash["Ances"] ):
        start,end = loc_iter.span()
        all_filter.extend( xrange(start,end ) )
    for loc_iter in re.finditer("\-+", cache_hash["mian"] ):
        start,end = loc_iter.span()
        all_filter.extend( xrange(start,end ) )    
    for data in [ cache_hash["Ances"],cache_hash["mian"]   ] :
        data = [ s for s in data  ]
        for i in all_filter:
            data[i]='_'
        data = ''.join(data)
        data = data.replace("_","")
        M_A.write(data+'\n')    
    #M_A.write(  cache_hash["Ances"]+'\n'+cache_hash["mian"]+'\n'   )
    M_A.close()
    mian_kaks = path+"/Mian_Anc.kaks"
    BASH.write("KaKs_Calculator  -i %s -o %s  -m NG \n"%(M_A.name, mian_kaks) )
    os.system(  "KaKs_Calculator  -i %s -o %s  -m NG 1>/dev/null 2>/dev/null "%(M_A.name, mian_kaks)  )
    #print("KaKs_Calculator  -i %s -o %s  "%(M_A.name, mian_kaks))
    if  os.path.getsize( mian_kaks  ):
        RAW = open( mian_kaks,'rU'  )
        RAW.next()
        #Ortholog_Pair[orthId]["KA/KS Harm"] = "Waiting"
        Ortholog_Pair[orthId]["KA/KS Harm"] = RAW.next().split("\t")[4]
    else:
        Ortholog_Pair[orthId]["KA/KS Harm"] = "Failed"
    
    yan_ances_name = path+"yan_vs_ances.axt"
    Y_A = open( yan_ances_name,'w' )
    Y_A.write("Yan_vs_Ancestor\n")
    all_filter = []
    for loc_iter in re.finditer("\-+", cache_hash["Ances"] ):
        start,end = loc_iter.span()
        all_filter.extend( xrange(start,end ) )
    for loc_iter in re.finditer("\-+", cache_hash["yan"] ):
        start,end = loc_iter.span()
        all_filter.extend( xrange(start,end ) )    
    for data in [ cache_hash["Ances"],cache_hash["yan"]   ] :
        data = [ s for s in data  ]
        for i in all_filter:
            data[i]='_'        
        data = ''.join(data)
        data = data.replace("_","")
        Y_A.write(data+'\n')

    Y_A.close()
    yan_kaks = path+"/Yan_Anc.kaks"
    os.system(  "KaKs_Calculator  -i %s -o %s -m NG 1>/dev/null 2>&1 "%(Y_A.name, yan_kaks)  ) 
    BASH.write("KaKs_Calculator  -i %s -o %s  -m NG\n"%(Y_A.name, yan_kaks) )
    if   os.path.getsize( yan_kaks ):
        RAW = open( yan_kaks,'rU'  )
        
        RAW.next()
        #Ortholog_Pair[orthId]["KA/KS Has"] = "Waiting!!"    
        Ortholog_Pair[orthId]["KA/KS Has"] = RAW.next().split("\t")[4]    
    else:
        Ortholog_Pair[orthId]["KA/KS Has"] ="Failed"

    



if __name__=="__main__":
    '# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
    usage='''usage: python %prog [options]

    Calcite KAKS '''
    parser = OptionParser(usage =usage )

    parser.add_option("-i", "--Table", action="store",
                      dest="Table",
                      type='string',
                      help="Ortholog Table")		




    parser.add_option("-o", "--Output", action="store", 
                      dest="output_Path",
                      help="Output Path prefix")
    (options, args) = parser.parse_args()
    outPATH = os.path.abspath(options.output_Path)+'/'
    if not os.path.exists(outPATH):
        os.makedirs(outPATH)
    BASH = open(outPATH+'run.bash','w')
    OrthoTable = pandas.read_table(options.Table)
    table_need = OrthoTable.loc[:,["Ortholog","H.armID","H.armSeq","H.asID","H.asSeq"]]
    Ortholog_Pair = Ddict()
    ALL_SEQ = open(outPATH+"Total_cds.fasta",'w')
    for i in xrange(0,len(OrthoTable)):
        table_data = OrthoTable.loc[i]
        
        orthId = table_data["OrthologID"]
        path_name = outPATH+orthId+'/'
        Ortholog_Pair[orthId][ "HarmId" ] = table_data[ "H.armID" ]
        Ortholog_Pair[orthId][ "HasId" ] = table_data[ "H.asID" ]
         
        if not os.path.exists(path_name):
            os.makedirs( path_name )
            
        RAW_SEQ = open(path_name+"/Unigene.fa",'w')
        RAW_SEQ.write(">"+table_data[ "H.armID" ]+'\n'+table_data[ "H.armSeq" ]+'\n'+ ">"+table_data[ "H.asID" ]+'\n'+table_data[ "H.asSeq" ]+'\n' )
        RAW_SEQ.close()
        ancestor = Pagan(RAW_SEQ.name)
        Ortholog_Pair[orthId][ "AncestorSeq" ] = ancestor
        
            
            
        
        if "Too " in ancestor:
            Ortholog_Pair[orthId][ "KA/KS Harm" ] = ancestor
            Ortholog_Pair[orthId][ "KA/KS Has" ] = ancestor     
        else:
            CDSPREDECTION  = open(path_name+"CDS_Predection.fa",'w'  )
            CDSPREDECTION.write( ">"+table_data[ "H.armID" ]+'\n'+table_data[ "H.armSeq" ]+'\n'+ ">"+table_data[ "H.asID" ]+'\n'+table_data[ "H.asSeq" ]+'\n'+">Ancestor_%s\n"%(orthId)+ancestor+'\n'   )
            CDSPREDECTION.close()
            cds_name = CdsFinder(CDSPREDECTION.name)
            if "Not " in cds_name:
                Ortholog_Pair[orthId][ "KA/KS Harm" ] = cds_name
                Ortholog_Pair[orthId][ "KA/KS Has" ] = cds_name
                
            else:
                KaksCal(cds_name)
                
    ALL_SEQ.close()
    KAKS_Result = open( outPATH+"/KAKS_Result.tsv",'w'   )
    KAKS_Result.write("OrthologID"+'\t')
    KAKS_Result.write(    '\t'.join(  sorted( Ortholog_Pair[ orthId ] )    ) +'\n'    )
    for orthId in Ortholog_Pair:
        KAKS_Result.write(orthId)
        for key in sorted( Ortholog_Pair[ orthId ]  ):
            KAKS_Result.write('\t'+Ortholog_Pair[ orthId ][ key ] )
        
        KAKS_Result.write('\n')