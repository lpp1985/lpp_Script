#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/3/30
"""

import unittest
from lpp import *
import os


if __name__ == '__main__':
    END = open("STATS.xls",'w')
    END.write("Sample\tReads\tMapped\tUniq Mapped\tMapped Perc%\tUnique Mapped Perc%\n")
    for a,b,c in os.walk("./"):
        for f in c:
            if f=="Star_OUTLog.final.out":
                RAW = open(a+'/'+f)
                data =RAW.read()
                all_name = glob.glob(a+'/*.bam')
                if all_name:
                    name = os.path.basename(all_name[0]).split(".")[0]
                    print(name)
                    all_reads = re.search("Number of input reads\s*\|\s+(\d+)",data).group(1)
                    uniq_mapped = re.search("Uniquely mapped reads number\s*\|\s+(\d+)",data).group(1)
		    multi_mapped = re.search("Number of reads mapped to multiple loci\s*\|\s+(\d+)",data).group(1)
		    all_mapped = int( uniq_mapped)+int(multi_mapped )
		    perc = float(all_mapped)*100/int(all_reads)

                    mapped_perc = 100*float( uniq_mapped ) /float( all_reads )
                    END.write("%s\t%s\t%s\t%s\t%.3f\t%.3f\n"%(name,all_reads,all_mapped,uniq_mapped,perc,mapped_perc))
                    
        
