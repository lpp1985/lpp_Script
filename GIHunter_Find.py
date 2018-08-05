#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/9
"""
from lpp import *

def Get_seq(data):
    seq_hash = {}
    RAW = fasta_check(  open(data,'rU')  )
    for t,s in RAW:
        name = t[1:].split()[0]
        seq_hash[name]=t.strip().split(" ",1)
        seq_hash[name].append(s)
    return seq_hash
if __name__ == '__main__':

    usage = "python2.7 %prog [options]"
    parser = OptionParser(usage =usage )
    parser.add_option("-i", "--InputPrefix", action="store",
                      dest="InputPrefix",

                      help="Contig")
    parser.add_option("-o", "--Output", action="store",
                      dest="OutputPrefix",

                      help="OutputPath")

    (options, args) = parser.parse_args()
    InputPrefix = options.InputPrefix
    Outprefix   = options.OutputPrefix
    # out_path = check_path( os.path.dirname(outprefix)  )


    ##Prepare Sequence
    GENOME = fasta_check(  open(InputPrefix+'.fna') )

    for t,s in GENOME:
        genomeseq  = re.sub("\s+", "", s)
    pep_hash = Get_seq(InputPrefix+'.faa')
    nucl_hash = Get_seq(InputPrefix+'.ffn')


    ##Run GIHunter##

    os.system(  "/home/nfs/SOFTWARE/Other/GIHunter/GIHunter  %(Input)s.fna  %(Input)s.ptt %(Input)s.rnt  %(Output)s "%(

        {
            "Input":InputPrefix,
            "Output":Outprefix,

        }
    )   
                )
    output_path = os.path.abspath(os.getcwd())+'/Genomic_Island/'+Outprefix+'/'

    if os.path.exists(output_path+'GIV/'):
        shutil.rmtree( output_path+'GIV/'  )
    #Prepare OUTPUT
    README = open( output_path+'Readme.txt','w' )
    README.write("""
使用GIHunter（http://www5.esu.edu/cpsc/bioinfo/software/GIHunter/）预测基因岛的结果
*_GIs.txt		GIHunter预测的原始结果（如果没发现基因岛的话，该文件为空）
*_GI.xls		基因岛的明细信息
Genome1_GI.stat	基因岛的长度统计信息
*_GIProtein.faa	基因岛内包含的蛋白质序列
*_GIGene.ffn		基因岛内包含的基因序列
*_GI.fa			基因岛完整序列
*_GI_Component.xls	基因岛包含基因的明细信息
    
    
    
    """)
    GISTAT =  open( output_path+Outprefix+'_GI.stat','w' )
    GISTAT.write( "\tName\tFrom\tEnd\tLength\n"    )
    GISEQ = open(output_path+Outprefix+'_GI.fa','w')
    GIGENE = open(output_path+Outprefix+'_GIGene.ffn','w')
    GIPROTEIN = open(output_path+Outprefix+'_GIProtein.faa','w')
    GENE_BELONG = open(output_path+Outprefix+'_GI_Component.xls','w')
    GENE_BELONG.write("Name\tBelongGI\n")
    GITSV = open(output_path+Outprefix+'_GI.xls','w')
    GITSV.write("\t".join( 
        [
            "Name","Kind","Function","Ref_Source" "Ref_Start","Ref_Stop" ,"Ref_Frame","Seq_Nucleotide",  "Seq_Nucl_Length" 
        ]  
        )+'\n'
            )
    gi_stat = []
    gi_location = {}

    #Parse GIHunter Result
    GIHUNTER_OUT = open(output_path+"out_GIs.txt" ,'rU' )
    gi_number = 0
    for line in GIHUNTER_OUT:
        gi_number+=1
        start,stop = line.strip().split()
        gi_name = Outprefix+"_GIsland%s"%( gi_number )
        gi_seq = genomeseq[int(start):int(stop)]
        GISEQ.write('>'+gi_name+'\n'+gi_seq)
        gi_length = int(stop)-int(start)+1
        gi_stat.append(gi_length)
        GISTAT.write('\t'+gi_name+'\t'+start+'\t'+stop+'\t'+str(gi_length)+'\n')
        gi_location[int(start)]=[ int(stop),gi_name   ]
        
        GITSV.write(  
        '\t'.join(
            [
                gi_name,
                "Genomic_Island",
                "Genomic_Island",
                Outprefix,
                start,
                stop,
                '+',
                gi_seq,
                str(gi_length)
            
            
            ]
        
        
            )+'\n'
        
    )

    PREDICTIONXLS = InputPrefix+'.xls'
    predict_frame = pd.read_table(PREDICTIONXLS)
    for i in xrange(0, len(predict_frame)):
        data =  predict_frame.loc[i]
        start = data["Ref_Start"]
        stop = data["Ref_Stop"]
        gi_belong = ""

        for gi_start,gi_list in gi_location.items():
            if start>= gi_start and start <= gi_list[0]:
                gi_belong = gi_list[1]
                break
            
            if stop>= gi_start and stop <= gi_list[0]:
                gi_belong = gi_list[1]
                break            
            
        if gi_belong:
            GENE_BELONG.write( data["Name"]+'\t'+gi_belong+'\n'  )
            
            GIGENE.write(
                
                nucl_hash[ data["Name"] ][0]+'___%s  '%(gi_belong)+\
                nucl_hash[ data["Name"] ][1]+'\n'+\
                nucl_hash[ data["Name"] ][2]
            
            
            
            
            )  
            if data["Name"] in pep_hash:
                GIPROTEIN.write(
                
                    pep_hash[ data["Name"] ][0]+'___%s  '%(gi_belong)+\
                    pep_hash[ data["Name"] ][1]+'\n'+\
                    pep_hash[ data["Name"] ][2]
                
                
                
                
                )               
    


    import numpy as np
    total_gi_length = sum(gi_stat)
    average_gi_length = numpy.average(gi_stat)
    GISTAT.write("AverageLength\t%s\n"%(average_gi_length))
    GISTAT.write("Total Length\t"+"("+"%s/%s) "%( total_gi_length,len(genomeseq) )+"%.2f" %(100.0* total_gi_length /len( genomeseq)  )+'%\n'   )
