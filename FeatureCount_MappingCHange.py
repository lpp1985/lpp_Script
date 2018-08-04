from lpp import *
mapping = File_dict( open(sys.argv[1])  ).read(1,2)
RAW = open(sys.argv[2],'rU')
END= open(sys.argv[3],'w')
for line in RAW:
	line_l = line.split("\t")
	if line_l[0] in mapping:
		line_l[0] = mapping[ line_l[0]  ]
	else:
		print(line_l[0])
	END.write('\t'.join(line_l)  )
