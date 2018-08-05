# usage python2.7 RAW_FILE   TARGET_FILE
from lpp import *
RAW = fasta_check(  open(sys.argv[1],'rU')  )
END = open(sys.argv[2],'w')
i=1
for t,s in RAW:
	s_new = re.sub('\s+',s)
	if len(s_new)<100:
		continue
	END.write(  '>%s\n'%(i)+s  )
	
	
	i=i+1