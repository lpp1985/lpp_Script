#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/9/3
"""
from lpp import *
import os
all_refseq = glob.glob("/home/nfs/SOFTWARE/Other/TransposonPSI_08222010/transposon_PSI_LIB/*refSeq")

RAW = fasta_check( open(sys.argv[1],'rU') )
TOP =open("Transpoase.tophit",'w')
ALL = open("Transpoase.allhit",'w')

TMP = os.path.abspath("./tmp/" )+'/'
def check_path(path):
	if not os.path.exists(path):
		os.makedirs(path)
	return os.path.abspath(path)+'/'
check_path(TMP)

if __name__ == '__main__':
	COMMAND = open(TMP+ "run.sh",'w')
	for t,s in RAW:
		name = t.split()[0][1:]
		name_path = check_path(TMP+name+'/')
		INPUTSEQ = open(name_path+name +".pep",'w')
		INPUTSEQ.write('>'+name+'\n'+s)
		INPUTSEQ.close()
		
		os.system("formatdb -i %s -p T"%( INPUTSEQ.name ))
		
		for e_ref in all_refseq:
			f_name = INPUTSEQ.name+"."+os.path.basename(e_ref).split(".")[0]
			COMMAND.write("blastpgp -i %s -d %s -R %s -j 1 -e 1e-5> %s && /home/nfs/SOFTWARE/Other/TransposonPSI_08222010/scripts/BPbtab <%s >%s\n"%(e_ref,INPUTSEQ.name, e_ref.rsplit(".")[0]+'.chk',f_name+'.psitblastn' ,f_name+'.psitblastn',f_name+'.btab' ) )
			
	os.system("cat %s|parallel -j 64"%(COMMAND.name) )
	for e_dir in glob.glob(TMP+'/*/'):
		all_bab = glob.glob(e_dir+'/*.btab')
		all_data = []
		for e_f in all_bab:
			RAW = open(e_f,'rU')
			for line in RAW:
				line_l = line.strip().split("\t")
				if len( line_l )<3:
					continue
				else:
					all_data.append(line_l)
		all_data = sorted(all_data,key=lambda x: float(x[-1]))
		if len(all_data)>0:
			TOP.write("\t".join(all_data[0])+'\n')
			for k in all_data:
				ALL.write("\t".join(k)+"\n")
		
		
