#!/usr/bin/python
from lpp import *
all_data = sys.argv[1:]
END = open("Maping.stast",'w')
END.write("Sample\tTotal reads\tMapped reads\tForward strand\tReverse strand\tFailed QC\tDuplicates\n")
for e_f in all_data:
	sample = os.path.basename(e_f).split(".")[0]
	RAW = open(e_f)
	RAW.next()
	RAW.next()
	RAW.next()
	RAW.next()
	print(RAW.next())
	END.write(sample)
	for line in RAW:
		if "Paired-end" in line:
			break
		END.write( "\t"+re.split("\:\s+",line.strip())[1].replace("\t",' ')  )
	END.write("\n")
