#!/usr/bin/python
import sys
RAW = open(  sys.argv[1],'rU')
LIST = open(sys.argv[2],'rU')
all_need = {}
for line in LIST:
	all_need[ line.split()[0]  ] = ""
print( all_need  )
END = open(sys.argv[3],'w')
END.write( RAW.next() )
for line in RAW:
	line_l = line.split("\t")
	if line_l[0]  in all_need:
		END.write(line)

