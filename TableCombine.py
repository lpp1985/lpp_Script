#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from lpp import *
from optparse import OptionParser
import os,string
def combine_xls( data1,data2,name1,name2   ):


    out_frame = pd.DataFrame.merge(data1, data2, left_on=name1,right_on = name2 , how='outer')
    return out_frame.drop_duplicates()



if __name__=="__main__":


    usage = '''usage: python2.7 %prog'''
    parser = OptionParser(usage =usage ) 
    
    parser.add_option("-n", "--Name", action="store", 
                      dest="Name", 

                      help="ColumName List (seperated by ';' )")
    parser.add_option("-o", "--outPath", action="store", 
                      dest="out", 
                      help="out")	



    (options, args) = parser.parse_args() 

    data_list = args
    name_list = options.Name.split('; ')
    out_frame = pd.read_table(data_list[0])
    for i in xrange(1,len(data_list)):
        out_frame = combine_xls( out_frame,pd.read_table(data_list[i]) , name_list[0],name_list[0]   )
    out_frame.to_csv(options.out,sep="\t",index=False)
    
    

    
