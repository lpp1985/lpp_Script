#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/6
"""
from optparse import OptionParser
import pandas as pd


dataframe = dataframe.drop("Nt_num",1  )
dataframe = dataframe.drop("Nt_num.1",1  )
dataframe = dataframe.drop("Nt_score",1  )
dataframe = dataframe.drop("Nt_query-frame",1  )  
dataframe = dataframe.drop("Nt_positive",1  )
dataframe = dataframe.drop("Nt_align-len",1  ) 
dataframe["Nt_Hit"] = dataframe["Nt_id"].astype("string")+" " + dataframe["Nt_def"].astype("string")
dataframe = dataframe.drop(["Nt_id","Nt_def","Nt_accession"],1  ) 
dataframe.rename( columns={"Nt_len":"Nt_Hit_len"},inplace=1     )
dataframe["Nt_Hit_Coverage"] = '('+(  abs(  dataframe["Nt_hit-from"]-dataframe["Nt_hit-to"]  )+1  ).astype("string")+'/'+dataframe["Nt_Hit_len"].astype("string")  +') '+ ( (  abs(  dataframe["Nt_hit-from"]-dataframe["Nt_hit-to"]  )+1  )/dataframe["Nt_Hit_len"] ).astype("string")

dataframe["Nt_Query_Coverage"] = '('+(abs(  dataframe["Nt_query-from"]-dataframe["Nt_query-to"]  )+1).astype("string")+'/'+dataframe["Nt_query-len"].astype("string")  +') '+ (100*(  abs(  dataframe["Nt_query-from"]-dataframe["Nt_query-to"]  )+1  )/dataframe["Nt_query-len"]).astype("string")

dataframe["Nt_identity"]  = 100*dataframe["Nt_identity"]/dataframe["Nt_query-len"]
dataframe2 = pd.DataFrame(dataframe,columns=['Name',u'Nt_Hit',"Nt_identity", 'Nt_query-len',u'Nt_Query_Coverage', u'Nt_query-from',u'Nt_query-to', u'Nt_Hit_len','Nt_Hit_Coverage',u'Nt_hit-frame',  u'Nt_hit-from',u'Nt_hit-to', u'Nt_evalue' , u'Nt_bit-score' ])



