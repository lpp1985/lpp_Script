#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from Dependcy import *
from optparse import OptionParser

if __name__=="__main__":
    config_hash = Config_Parse()

    usage = '''usage: python2.7 %prog'''
    parser = OptionParser(usage =usage ) 
    parser.add_option("-p", "--PEP", action="store", 
                      dest="PEP", 
                      default = "",
                      help="protein file")
    parser.add_option("-n", "--NUL", action="store", 
                      dest="NUL", 
                      help="necleotide file")		

    parser.add_option("-o", "--end", action="store", 
                      dest="output_prefix", 
                      help="output_prefix")

    parser.add_option("-e", "--evalue", action="store", 
                      dest="evalue", 
                      help="evalue cutoff")
    parser.add_option("-k", "--ko", action="store_true", 
                      dest="KO", 
                      default=False,
                      help="evalue cutoff")    



    (options, args) = parser.parse_args() 
    proteinseq = options.PEP
    if not proteinseq:
        proteinseq = options.NUL
    if not options.NUL:
	options.NUL=proteinseq    
    FASTA = fasta_check(  open(proteinseq,'rU')  )
    sequence = FASTA.next()[-1]
    blast_type = Nul_or_Protein(sequence)
    output_prefix = os.path.abspath(  options.output_prefix )
    
    out_put_path = os.path.split(output_prefix)[0]+'/'

    tag = "%s"%( os.getpid() )
    if not os.path.exists( out_put_path ):
        os.makedirs( out_put_path )
    end_list = glob.glob(out_put_path+'/*.zip')
    if end_list:
        sys.exit()
    README = open(out_put_path+"/Readme.txt",'w')
    README.write(
"""
该文件夹放置KEGG通路分析结果，我们将所得的序列比对到KEGG数据库并进一步映射到KO（KEGG Orthology），并通过KO映射到Pathway。附件说明如下：
*_pathway.tsv\t基因映射到KEGG KO 和Pathway的明细，用Excel打开
*_AlignKEGG.tsv\t基因序列与KEGG数据库比对结果，用excel打开
*_PathwayCategory_Stats.stat\tKEGG每一个功能模块的Pathway映射到的基因个数统计,用excel打开
stat.*\tPathway分析可视化结果，提供tiff和PDF两个版本
*.gz\tPathway分析结果的可视化结果，请解压，有两个文件夹，其中doctree文件夹位搜索索引，无需打开。Pathway文件夹下是一个网站，请点击index.html观看，每一个Pathway如果有基因被比对上，后面会出现all字样。
*.R\t可视化画图脚本，用R语言运行，您可以根据需要自行调整。



"""
    
    
    
    
    )        
    diamond_result = output_prefix+'_AlignKEGG.tsv'
    
        
    if options.KO:
        error = RunDiamond(proteinseq,options.evalue, blast_type,"ko",diamond_result)
    else:
        
        error = RunDiamond(proteinseq,options.evalue, blast_type,"kegg",diamond_result)
    #error=""
    if error:
        print( colored("%s 's KEGG process in Diamond of kegg is error!!","red") )
        print(colored( error,"blue"  ))
        print(  "##############################################"   )

        sys.exit()

 
    blast_mapping_command = config_hash["Utils"]["gapmap"]+'/blast_sql.py -f %(diamond)s   -r  %(diamond)s   -1 forward -2 forward  -n Forward -N Reverse -p %(pep)s -d %(dna)s -x %(tag)s -q'%(
        {
            "diamond":diamond_result,
            "pep":proteinseq,
            "dna":options.NUL,
            "tag":tag
        }
    )
    print( blast_mapping_command  )
    blast_mapping_process = subprocess.Popen( blast_mapping_command.split(),stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
    blast_mapping_process.communicate()




    source_location = out_put_path+'/source/'

    source_command = config_hash["Utils"]["gapmap"]+"/Mapping_sql.py -o %s -d %s"%(
        source_location,
        tag,
        ) +" && "+config_hash["Utils"]["gapmap"]+"/Show_all_Mapping.py -o %s  -d %s"%(
            output_prefix,
            tag,
        )
    print(source_command  )
    os.system(  source_command )

    pathway_detail_frame = pd.read_table( "%s_detail.tsv"%(  output_prefix  )   )
    pathway_detail_frame = pd.DataFrame(pathway_detail_frame,columns=pathway_detail_frame.columns[:-1])
    os.remove( "%s_detail.tsv"%(  output_prefix  ) )
    
    pathway_stats_command = config_hash["Utils"]["gapmap"]+"/Pathway_stats.py %(name)s %(out)s_PathwayCategoery2.tsv %(out)s_PathwayCategory_Stats.stat  %(out)s_PathwayCategoery.tsv"%(
                {
                    "name":tag,
                    "out":output_prefix
                }
                )
    pathway_stats_process = subprocess.Popen( pathway_stats_command.split(),stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
    pathway_stats_process.communicate()
    pathway_category_frame = pd.read_table("%s_PathwayCategoery2.tsv"%(output_prefix) )
    
    pathway_result_frame = pd.merge( pathway_detail_frame,pathway_category_frame,left_on='Name', right_on='Name', how='outer' )
    pathway_result_frame.to_csv( "%s_pathway.tsv"%(output_prefix),sep="\t",index=False  )
    os.remove("%(out)s_PathwayCategoery2.tsv"%(
    {
        "out":output_prefix
    }
    )  
              )


    pathwaydraw_command = "Pathway_Draw.py   -i %s_pathway.tsv  -o %s -r %s"%(
        output_prefix,
        out_put_path+"stats",
        out_put_path+'Draw.R',
    )
    pathwaydraw_process = subprocess.Popen(  pathwaydraw_command.split(),stderr= subprocess.PIPE,stdout=  subprocess.PIPE  )
    stdout,stderr = pathwaydraw_process.communicate()	

    pathway_category_frame = pd.read_table("%s_PathwayCategoery.tsv"%(output_prefix) )

    pathway_result_frame = pd.merge( pathway_detail_frame,pathway_category_frame,left_on='Name', right_on='Name', how='outer' )
    pathway_result_frame.to_csv( "%s_pathway.tsv"%(output_prefix),sep="\t",index=False  )
    
    
    database_generate_commandline = " KEGG_Database.py  -a %s  -p %s  -d %s_KEGG.xls  "%(  diamond_result,"%s_pathway.tsv"%(output_prefix) ,output_prefix   )
    os.system(database_generate_commandline)
    
    
    
    make_commandline = config_hash["Tools"]["sphinx"] +" -j 32 -Q -b html -d %(out)s/doctrees %(location)s %(out)s/pathway && zip  %(out)s.zip %(out)s -rq && rm %(out)s -rf &&rm %(location)s -rf"%(
        {
            "out":out_put_path+"Pathway",
            "location":source_location
        }
    ) 
    print( make_commandline ) 
  
    os.system(make_commandline)
    os.remove( config_hash["Utils"]["gapmap"]+'/'+tag )






