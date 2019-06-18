#!/usr/bin/python
from lpp import *
all_data = glob.glob("*.stats")
END = open( "Stats.xls",'w'  )
i=0
for f in all_data:
	RAW = open(f,'rU'  )
	if i==0:
		i=1
		END.write("Sample\t"+ RAW.next())
		
	else:
		RAW.next()
	name = f.rsplit(".",1)[0]
	END.write( name+'\t'+RAW.next() )

