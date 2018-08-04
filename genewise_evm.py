#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/7/21
"""
from lpp import *
RAW = open(sys.argv[1],'rU')
END = open(sys.argv[2],'w')
all_has = {}
OLD_ID = ""
i=0
cache = []
for line in RAW:
    
    
    ID = re.search("ID=([^\;]+)",line).group(1)
    
    line_l = line.strip().split("\t")
    
    if not OLD_ID:
        
        OLD_ID= ID 
       
        
    if ID != OLD_ID :
        cache = sorted(  cache,key = lambda x: int(x.split("\t")[3]  )  )
        old_data = cache[0].split("\t")[:-1]
        attribute = "ID=exons.gene%s;Parent=Protein%s\n"%(i,i)
        end = cache[-1].split("\t")[4]
        old_data[4]=end

    
        old_data[2] = "gene"
        END.write("\t".join(old_data))
        attribute = "ID=gene%s;Name=gene%s\n"%(i,i)
        END.write('\t'+attribute)
        old_data[2] = "mRNA"
        attribute = "ID=Protein%s;Parent=gene%s\n"%(i,i)
        END.write("\t".join(old_data)+'\t'+attribute)
        i+=1
        END.write(  "".join(cache) )
        
        OLD_ID=ID
        
        
        cache = []
        attribute = "ID=exons.gene%s;Parent=Protein%s\n"%(i,i)
        line_l[2]="exon"
        cache.append(   '\t'.join(line_l[:-1])+'\t'+attribute  )
        attribute = "ID=cds.gene%s;Parent=Protein%s\n"%(i,i)
        line_l[2]="cds"
        cache.append(   '\t'.join(line_l[:-1])+'\t'+attribute  )
    else:
        
        attribute = "ID=exons.gene%s;Parent=Protein%s\n"%(i,i)
        line_l[2]="exon"
        
        cache.append(   '\t'.join(line_l[:-1])+'\t'+attribute  )
        attribute = "ID=cds.gene%s;Parent=Protein%s\n"%(i,i)
        line_l[2]="cds"   
        cache.append(   '\t'.join(line_l[:-1])+'\t'+attribute  )

       
    


    
    
        
        
        
else:
    cache = sorted(  cache,key = lambda x: int(x.split("\t")[3]  )  )
    old_data = cache[0].split("\t")[:-1]
    attribute = "ID=exons.gene%s;Parent=Protein%s\n"%(i,i)
    end = cache[-1].split("\t")[4]
    old_data[4]=end


    old_data[2] = "gene"
    END.write("\t".join(old_data))
    attribute = "ID=gene%s;Name=gene%s\n"%(i,i)
    END.write('\t'+attribute)
    old_data[2] = "mRNA"
    attribute = "ID=Protein%s;Parent=gene%s\n"%(i,i)
    END.write("\t".join(old_data)+'\t'+attribute)
    i+=1
    END.write(  "".join(cache) )

    
        