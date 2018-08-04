#!/usr/bin/python
from lpp import *
RAW = open(sys.argv[1],'rU')
data = re.split("\n\n+", RAW.read() )
END = open("emv_out.gff3",'w')
for e_b in data:
	i=0
	for line in e_b.split("\n"):
		line = re.sub("(?=\S+)\.",'_',line)
		line_l = line.split("\t")
		if line_l[2] =="CDS":
			i+=1
			line_l[-1] = line_l[-1].replace(";",".CDS%s;"%(i))
			
		END.write("\t".join(line_l)+'\n')
	END.write("\n\n")
