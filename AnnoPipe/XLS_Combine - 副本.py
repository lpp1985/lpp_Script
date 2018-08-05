#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from Dependcy import *
from optparse import OptionParser
import os,string
def combine_xls( data_list   ):
    out_frame = pd.read_table(data_list[0]).drop_duplicates()

    for each_data in data_list[1:]:

        new_frame = pd.read_table(each_data).drop_duplicates()
        on_need = list(out_frame.columns  & new_frame.columns)

        out_frame = pd.DataFrame.merge(out_frame, new_frame, on=on_need, how='outer')
    return out_frame.drop_duplicates()
def Rawcombine_xls( data_list   ):
    out_frame = data_list[0].drop_duplicates()

    for new_frame in data_list[1:]:

        on_need = list(out_frame.columns  & new_frame.columns)

        out_frame = pd.DataFrame.merge(out_frame, new_frame, on=on_need, how='outer')
    return out_frame.drop_duplicates()

if __name__=="__main__":
    config_hash = Config_Parse()

    usage = '''usage: python2.7 %prog'''
    parser = OptionParser(usage =usage ) 
    
    parser.add_option("-i", "--inputpath", action="store", 
                      dest="inputpath", 
                      default = "./",
                      help="input path")
    parser.add_option("-o", "--outPath", action="store", 
                      dest="outputpath", 
                      help="outpath")	
    parser.add_option("-g", "--GFFANNO", action="store", 
                      dest="gff", 
                      default="",
                      help="Annotation from gff file or Raw Unigene Sequence!")	    



    (options, args) = parser.parse_args() 

    out_put_path = os.path.abspath(  options.outputpath )+'/'

    if not os.path.exists( out_put_path ):
        os.makedirs( out_put_path )
    README = open(out_put_path+"/Readme.txt",'w')
    README.write(
"""
该文件夹放置所有注释分析的表格，分成两类，第一类是按照不同的数据库进行分类，每一个子文件夹放置的excel文件包含该数据库下所有样本的注释结果。第二类是按照样本来源分类，每一个子文件夹下放置不同染色体的基因注释信息附件说明如下：
Database文件夹\t不同数据库注释的汇总文件
Choromsome文件夹\t按照不同染色体进行的分类汇总文件
All_HasAnnotation.xlsx		所有能够注释的基因的注释明细表
GeneFeature+Annotation.xlsx	注释的基因信息和基因序列等信息的总表




""")
    category_hash = Ddict()
    chrosome_hash = Ddict()
    for base_path,dir_name,file_list in os.walk(options.inputpath):
        for e_f in file_list:
            if e_f.endswith(".xls"):
                chrosmome,category=  base_path.rsplit('/',2)[-2:]
                category_hash[category][chrosmome] = base_path+'/'+e_f
                chrosome_hash[chrosmome][category] = base_path+'/'+e_f
                
    chrosome_dir = out_put_path+"/Choromsome/"
    category_dir = out_put_path+"/Database/"
    for e_path in [ chrosome_dir, category_dir  ]:
        if not os.path.exists(e_path):
            os.makedirs(e_path)
            
    category_Excel = pd.ExcelWriter(category_dir+'CategoryAnnotation.xlsx', engine='xlsxwriter')
    STAT = open(category_dir+'/stat.tsv','w')
    STAT.write("Database\tHitGeneNumber\tTotalGeneNumber\tPerc(%)\n")
    total_excel = []
    database_data  = {}
    stat_result = {}

    for category in category_hash:
        all_excel = []
        
        for chrosome in category_hash[category]:
            all_excel.append(  category_hash[category][chrosome]  )
        total_excel.extend(all_excel)
        result_frame = combine_xls(all_excel)
        stat_result[category] = [category,len(result_frame["Name"])  ]
        # STAT.write(category+'\t%s\t%.2f\n'%(len(result_frame["Name"] ) ,100.0* len(result_frame["Name"] )  /   ) )
        
        result_frame["from"] = result_frame["Name"].str.split('_',1).str.get(0)
        
        result_frame["id"] = result_frame["Name"].str.split('_',1).str.get(1)
        result_frame =result_frame.sort(["from",'id'],axis=0)
        result_frame = result_frame.drop(["from",'id'],axis=1)
        
        result_frame.to_excel( category_Excel,category ,index=False   )
        if category!="Nt":
            database_data[category] = result_frame["Name"]
    category_Excel.save()
    VENN_R = open( category_dir+"/Draw.R",'w'   )
    VENN_R.write(
        """
require("VennDiagram")
temp = venn.diagram(
     x = list(
 
        """)
    
    
    
    end_list = []
    for category ,data in database_data.items():
        end_list.append(  """    %s=c(%s)
        
"""%(
    category,
    ','.join(["'"+x+"'"  for x in data])
    
   
   )
   )
    VENN_R.write(",".join(end_list))

    VENN_R.write("""),
    filename = NULL,
	col = "black",
	lty = "solid",
	lwd = 4,
	fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
	alpha = 0.50,
	label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
	cex = 2.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
	cat.cex = 2.5,
	cat.fontfamily = "serif"
	)    
pdf("%s")
grid.draw(temp)    
dev.off()    
    """%(
           category_dir+'/stat.pdf'
           
       
       )
       
       
       )
    VENN_R.close()
    os.system("Rscript %s "%( VENN_R.name  ))
    
    
    chrosome_Excel = pd.ExcelWriter(chrosome_dir+'ChorosomeAnnotation.xlsx', engine='xlsxwriter')
    all_result = []
    for chrosome in chrosome_hash:
        all_excel = []
        for category in chrosome_hash[chrosome]:
            all_excel.append(  category_hash[category][chrosome]  )
        
        result_frame = combine_xls(all_excel)
        result_frame["from"] = result_frame["Name"].str.split('_',1).str.get(0)
    
        result_frame["id"] = result_frame["Name"].str.split('_',1).str.get(1)
        result_frame =result_frame.sort(["from",'id'],axis=0)
        result_frame = result_frame.drop(["from",'id'],axis=1)        
        result_frame.to_excel( chrosome_Excel,chrosome,index=False    )
        all_result.append(result_frame)
    chrosome_Excel.save()
    
    all_resultframe = Rawcombine_xls(all_result)

    
    all_resultframe["from"] = all_resultframe["Name"].str.split('_',1).str.get(0)
    all_resultframe["id"] = all_resultframe["Name"].str.split('_',1).str.get(1)
    all_resultframe =all_resultframe.sort(["from",'id'],axis=0)
    all_resultframe = all_resultframe.drop(["from",'id'],axis=1)
    all_resultframe = all_resultframe.drop_duplicates()
    all_resultframe.to_csv( out_put_path+"All_HasAnnotation.xlsx" ,index=False ,sep='\t'  )
    stat_result["Total"] = [ "Total",len( all_resultframe["Name"]  ) ]
    if options.gff:
        DATA = open(options.gff ,'rU')
        line = DATA.next()
        if line.startswith('>'):
            total_number = 0
            SEQ = fasta_check(open(options.gff ,'rU') )
            for t,s in SEQ:
                total_number +=1
            
        else:
            all_gff_frame= pd.read_table( options.gff  )
            total_resultframe = pd.DataFrame.merge(all_gff_frame, all_resultframe, left_on='Name', right_on='Name', how='outer')
            total_number = len( total_resultframe["Name"]   )
            
            total_resultframe.to_csv(out_put_path+"GeneFeature+Annotation.xlsx",index=False,sep="\t"   )
        for key in stat_result: 
            if key !="Total":
                category,genenumber = stat_result[key]
                perc = 100.0*genenumber/total_number
                STAT.write(  "\t".join( [ category, str(genenumber),str(total_number) ,"%.2f"%(perc)     ] ) +'\n'   )
        
        category,genenumber = stat_result["Total"]
        perc = 100.0*genenumber/total_number
        STAT.write(  "\t".join( [ category, str(genenumber),str(total_number) ,"%.2f"%(perc)     ] ) +'\n'   )        
        
    
    

    