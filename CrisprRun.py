#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/8
"""
from lpp import *

import pandas as pd
if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] 
		     Pilr-CR Crispr Finding !!!
		     '''
    
    parser = OptionParser(usage =usage )    
    parser.add_option("-g", "--Genome", action="store",
                      dest="Genome",
                      help="Genome Sequence")


    parser.add_option("-o", "--OUTPUT", action="store",
                      dest="outputprefix",
                      help="Output Path") 

    (options, args) = parser.parse_args()
    outputprefix = options.outputprefix
    output_path = check_path(os.path.dirname(outputprefix) )    
    END = open(output_path+'/Readme.txt','w')
    END.write(
"""
使用Pile-CR预测
  *_DP.fa Crispr串联重复序列
  *_Spacer.fa   Crispr中Spacer序列
  *_Spacer.xls   Crispr中Spacer详细注释结果和基因组位置
  *_NtAlign.tsv   Crispr中Spacer的序列注释结果
  *_RAW.txt     Pilercr原始分析结果
  


"""    
    
    
    
    )
    
    Genome = os.path.abspath(options.Genome)
    TMP_INPUT = open( output_path+"%s.contigs"%(os.getpid()),'w' )
    seq_hash = {}
    for t,s in fasta_check( open(Genome,'rU') ):
        t = re.sub( "_+$","", t.strip().split("|")[-1])
        if t.startswith('>'):
            t = t[1:]
            
        s1 = re.sub("\s+", "", s)
        seq_hash[t]=s1
        TMP_INPUT.write('>'+t+'\n'+s)
        
    TMP_INPUT.close()
        
    
    tmp_name = "%s_RAW.txt"%(outputprefix)
    os.system ("pilercr -in %s  -out %s -seq %s_DP.fa -quiet -noinfo -trimseqs"%(TMP_INPUT.name,tmp_name,outputprefix)) 

    os.remove( TMP_INPUT.name)

    RAW = open(tmp_name,'rU')
    total_data = RAW.read()
  
    
    ###Justify if has Crispr Sequence!!
    if "0 putative CRISPR arrays found." in total_data:
        END = open(output_path+'Result.txt','w')
        END.write("0 putative CRISPR arrays found.")
        os.remove("%s_DP.fa"%( outputprefix  ))
        sys.exit()
        
    

    #####Rename DP Sequence ######
    DP_SEQ =    fasta_check( open( "%s_DP.fa"%(outputprefix) ,'rU' )  )
    Crispr_id = 0
    DP_NEW =  open( "%s_DP2.fa"%(outputprefix) ,'w' )
    for t,s in DP_SEQ:
        Crispr_id +=1
        genome_name = t.strip().split("[")[0]
        DP_NEW.write(genome_name+'_CrisprEmement%s\n'%(Crispr_id)+s)
    
    shutil.move(  DP_NEW.name,"%s_DP.fa"%(outputprefix))    
    
    
    
    
    ####################SPACER PARSE#########################################
    data_block = re.split( "\n(?:DETAIL REPORT|SUMMARY BY [A-Z]+)\n",total_data    )
    align_detail = data_block[1]
    align_detail = align_detail.strip()
    detail_list = re.split("\n{3}",align_detail)
    SPACER_SEQ = open(outputprefix+'_Spacer.fa','w')
    SPACER_TSV = open(outputprefix+'_Spacer.xls','w')
    SPACER_TSV.write( 
        '\t'.join( 
            ["Name","Kind","Function","Ref_Source","Ref_Start","Ref_Stop","Ref_Frame"	,"Seq_Nucleotide",	"Seq_Nucl_Length"]
            )+'\n' 
                      )
    AlignList = namedtuple("Align","Pos,RepeatLength,iden,SpacerLength,Left,Repeat,Spacer")
    for each_detail in detail_list:
        crispr_number = re.search("Array\s+(\d+)",each_detail).group(1)
        seq_name = re.search(">(\S+)",each_detail).group(1)
        data_line_list = each_detail.split("\n")
        data_line_list = data_line_list[5:]
        spacerid = 0 
        for key in data_line_list:
            align_list = key.split()
            if key.startswith("=="):
                break
            if len(align_list)!=7:
                continue
            spacerid+=1
            align_list = AlignList._make(align_list)
            startpos = int(align_list.Pos)
            repeat_length = len(align_list.Repeat)
            spacer_start = startpos+repeat_length
            spacer_end = spacer_start+int(align_list.SpacerLength)
            spacer_name = seq_name+"_Crispr%sSpacer%s"%(crispr_number,spacerid)
            
            SPACER_TSV.write(
                "\t".join(
                    [
                        spacer_name,
                        "CrisprSpacer",
                        "CrisprSpacer",
                        seq_name,
                        str(spacer_start),
                        str(spacer_end),
                        "+",
                        align_list.Spacer,
                        
                        str(align_list.SpacerLength),
   
                    ]            
                )+'\n'
            )
            SPACER_SEQ.write('>'+spacer_name+'\n'+align_list.Spacer+'\n')
            
            
            
            
        
            
    



    ################Element Parse #########################
    element_block =data_block[-1].strip()
    element_list = element_block.split("\n")[4:]
    ELEMENT_SEQ = open(  "%s_CrisprElement.fa"%( outputprefix  ),'w')
    ELEMENT_TSV = open(  "%s_CrisprElement.xls"%( outputprefix  ),'w')
    ELEMENT_STAT = open(  "%s_CrisprElement.tsv"%( outputprefix  ),'w')
    title_data = element_block[0]
    title_list= ["Name","Ref_Source" ,"Ref_Start","Ref_Stop" ,"Ref_Frame","Seq_Nucleotide",  "Seq_Nucl_Length"  , "Copies" , "Repeatlength" , "Spacerlength" ,"Distance" , "DP_Consensus"]
    ELEMENT_TSV.write( 
        "\t".join( 
            [
                "Name","Kind","Function","Ref_Source" ,"Ref_Start","Ref_Stop" ,"Ref_Frame","Seq_Nucleotide",  "Seq_Nucl_Length" 
            ]  
        )+'\n'
                       )
    ELEMENT_STAT.write( "\t".join(  title_list )+'\n')
    ElememtList = namedtuple(
    "Length",
    "Number,Source,Start,Length,Copies,DpRepeatLength,SpacerLength,Distance,DPConsensus")
    
    for each_element in element_list:
        data_list = each_element.split()
        if len(data_list)!=9:
            data_list.insert(-1,"")
        element_detail_list = ElememtList._make( data_list )
        elementname = element_detail_list.Source+"_CrisprElement"+element_detail_list.Number
        elementstop = str( int(element_detail_list.Start)+int(element_detail_list.Length) )
        elment_seq = seq_hash[ element_detail_list.Source  ][  int(element_detail_list.Start)   :int(elementstop) ]
      
        ELEMENT_SEQ.write('>'+elementname+'\n'+elment_seq+'\n')
        ELEMENT_TSV.write(
            '\t'.join(
                [
                    elementname,
                    "CrisprElement",
                    "CrisprElement",
                    element_detail_list.Source,
                    element_detail_list.Start,
                    elementstop,
                    '+',
                    elment_seq,
                    element_detail_list.Length
                    
                
                
                ]    
            )+'\n'        
        )
        
        ["Name","Ref_Source" "Ref_Start","Ref_Stop" ,"Ref_Frame","Seq_Nucleotide",  "Seq_Nucl_Length"  , "Copies" , "Repeatlength" , "Spacerlength" ,"Distance" , "DP_Consensus"]
        ELEMENT_STAT.write(
            '\t'.join(
                [
                    elementname,
                    element_detail_list.Source,
                    element_detail_list.Start,
                    elementstop,
                    '+',
                    elment_seq,
                    element_detail_list.Length,
                    element_detail_list.Copies,
                    element_detail_list.DpRepeatLength,
                    element_detail_list.SpacerLength,
                    element_detail_list.Distance,
                    element_detail_list.DPConsensus
            
            
            
                ]    
                )+'\n'             
            
        
        
        
        )
        
    
    ## SPACE NT ALignment ##
    os.system( "Nt_Annotation.py  -i  %s -o  %s_NtAlign.tsv  -e 1e-5 "%(SPACER_SEQ.name,outputprefix)  )
    space_raw_frame = pd.read_table(SPACER_TSV.name)
    if os.path.getsize(  "%s_NtAlign.tsv" %(  outputprefix )   ):
        space_anno_frame = pd.read_table(   "%s_NtAlign.tsv" %(  outputprefix )    )
        space_new_frame = pd.DataFrame.merge(  space_raw_frame, space_anno_frame,left_on="Name",right_on="Name",how='outer'  )
        space_new_frame.to_csv( SPACER_TSV.name,sep="\t",index=False    )
    else:
        space_raw_frame.to_csv( SPACER_TSV.name,sep="\t",index=False )