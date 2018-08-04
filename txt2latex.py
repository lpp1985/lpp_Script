#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/5/9
"""
import re
from optparse import OptionParser
def trans(data):
	data_list = data.rstrip().split('\t')
	

	for i in xrange( 0 , len(data_list) ):
		if re.match(  "^\d+\.\d+$"  ,  data_list[i]   ):
			data_list[i] = "%.2f"%(  float(data_list[i])  )
			
	data = "\t".join(data_list)+'\n'
	data=data .replace( "\\","\\\\")
	data = re.sub(  "(?=[^\\\])%","\%",data )
	data = re.sub(  "(?=[^\\\])_","\_",data )
	data=data.replace("\t","&")
	data = data.replace("\n","\\\\")
	return data
	
if __name__ == "__main__":
	usage = '''usage: python2.7 %prog [options] 
'''
	parser = OptionParser(usage =usage )

	parser.add_option("-i", "--Input", action="store",
	                  dest="inputData",
	                  help="Input Data")	
	parser.add_option("-o", "--Output", action="store",
	                  dest="Output",
	                  help=" Output")		
	parser.add_option("-c", "--Caption", action="store",
	                  dest="Caption",
	                  help=" Caption")			
	(options, args) = parser.parse_args()
	inputData = options.inputData
	Output = options.Output
	Caption = options.Caption
	OUTPUT = open(  Output,'w' )
	
	RAW = open(inputData,'rU')
	title =RAW.next()
	title_l =  title.split("\t")
	format_tag =["c"]*len(  title_l) 
	for i in xrange(0,len(title_l)):
		if len(title_l[i])>15:
			format_tag[i]="X"
	
	for line in RAW:
		line_l = line.strip().split("\t")
		for i in xrange(0,len(line_l)):
			if len(line_l[i])>200:
				format_tag[i]="X"
				
	if len(format_tag)>8:
		OUTPUT.write("""\\begin{sidewaystable}\n
		\\begin{center}
	\\small
	\\setlength\\tabcolsep{1pt}\n""")
	else:
		OUTPUT.write("""\\begin{table}[H]\n\\begin{center}\n""")
	format_tag="{"+"".join(format_tag)+'}'
	
	OUTPUT.write("""
    \\caption{%s}
    
       \\begin{tabularx}{\\textwidth}"""%(Caption))	
	
	
	OUTPUT.write	(format_tag+'\n')
	RAW = open(inputData,'rU')
	title = trans(RAW.next())
	OUTPUT.write("\\hline\n"+title+'\n \\hline\n')
	for line in RAW:
		line = trans(line)
		OUTPUT.write(line+'\n')
	OUTPUT.write("\\hline\n\\end{tabularx}\n\\end{center}\n")
	if len(format_tag)>8:
		OUTPUT.write("\\end{sidewaystable}\n")
	else:
		OUTPUT.write("\\end{table}\n")
		
	
