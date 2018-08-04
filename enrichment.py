#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/4
"""

from numpy.random import hypergeometric
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects import r
from math import factorial

stats = importr('stats')
qval = importr("qvalue")

def enrichment_analysis(diff_gene,total_anno_gene,sample_size,sampled_diff):
    #超几何分布检验该通路是否显著变化




    p_value = r('1- phyper(%s, %s, %s, %s)'%(sampled_diff-1,diff_gene,total_anno_gene-diff_gene,sample_size))[0]

    return p_value
def fdr(p_value_list):    
    p_adjust = stats.p_adjust(FloatVector(p_value_list), method = 'BH')
    q_adjust = qval.qvalue(p_adjust)
    p_adjust_value = [k for i,k in p_adjust.items()]
    for i in q_adjust.items():
        if i[0]=="qvalues":
            q_adjust_value = i[1]
            
  
    return p_adjust_value,q_adjust_value
