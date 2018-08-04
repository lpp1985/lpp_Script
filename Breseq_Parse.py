#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/6/30
"""
from bs4 import BeautifulSoup
from lpp import *



usage = '''usage: python2.7 %prog'''
parser = OptionParser(usage =usage ) 



parser.add_option("-i", "--Html", action="store", 
                      dest="Index", 
                      help="Index.html")



parser.add_option("-k", "--SnpResult", action="store", 
                      dest="SnpResult", 
                      help="Snp Result")


parser.add_option("-m", "--IndelResult", action="store", 
                      dest="IndelResult", 
                      help="Indel Result")

(options, args) = parser.parse_args() 



html = open( options.Index, 'rU').read().replace("&#8209;", "-")
soup = BeautifulSoup(html)
all_table = soup.find_all("table", attrs={"cellspacing": "1"})
end_table = ""
result_table = []
for each_table in all_table:
	end_table = ""
	for each_tr in each_table.children:
		if each_tr.name != 'tr':
			continue
		end_list = []
		for each_th in each_tr.children:
			if each_th.name in[ "th", "td"]:
				end_list.append(each_th.text.strip())
		end_table += '\t'.join(end_list) + '\n'
	result_table.append(end_table)
		
END = open( options.SnpResult, 'w')
END.write(result_table[0].encode("utf-8"))
END = open(options.IndelResult, 'w')
END.write(result_table[1].encode("utf-8"))
