#!/usr/bin/env python
from lpp import *
for a,b,c in os.walk("./"):
	for f in c:
		if ".gz" in f:
			#print(a+"/"+re.sub("_FDMS\S+_","_", f.replace("smbh-","")    ))
			os.rename( a+"/"+f ,a+"/"+re.sub("_FDMS\S+_","_", f.replace("smbh","")  )  )
