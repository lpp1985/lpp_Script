#!/usr/bin/python
from XLS_Combine import combine_xls
from lpp import *
data_list = sys.argv[1:]
END = "Combined.xls"
out_frame = combine_xls( data_list )
out_frame.to_csv(  END,index=False,sep="\t")
