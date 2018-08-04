#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/12/15
"""
from lpp import *
from optparse import OptionParser
usage = '''usage: python2.7 %prog [options] 
         parse eggNOG data

        '''
parser = OptionParser(usage =usage )    
parser.add_option("-i", "--Input", action="store",
                dest="name",
                help="Input sequence name")
parser.add_option("-o", "--out", action="store",
                dest="outprefix",
                help="output path")


if __name__ == '__main__':
	(options, args) = parser.parse_args()
	seqname = options.name
	prefix = os.path.basename(seqname)
	outpath = options.outprefix+'/'+prefix
	outpath = check_path(outpath)
	print( " phage_finder_v2.0.sh  %s  %s"%(seqname,outpath) )
	outputname = outpath+prefix+'_Phage'
	os.system( " phage_finder_v2.0.sh  %s  %s"%(seqname,outpath)  )
	README = open(outpath+'/Readme','w')
	README.write("""
使用PhageFinder进行前噬菌体寻找。结果如下：
*.con 前噬菌体序列
*.xls前噬菌体序列的详细详细信息表格
*.pep前噬菌体内包含的蛋白
*.seq 前噬菌体内包含的基因
	
	
	""")
	for e_f in glob.glob(outpath+"*.*"):
		if e_f.endswith(".hmm3") or e_f.endswith(".out") or e_f.endswith(".log"):
			os.remove(e_f)
		else:
			path,file_name = os.path.split(e_f)
			appendix = file_name.rsplit(".",1)[-1]
			shutil.move(e_f,outputname+'.'+appendix)
			
	if os.path.getsize(outputname+'.pep'):
		
		XLS = open(outputname+"_PhageGene.xls",'w')
		XLS.write("Name\tBelongToPhage\n")
		TMP = open(outputname+'.tmp','w')
		for t,s in fasta_check(  open(outputname+'.seq'   )   ):

			seq_name,annotation = t[1:].strip().split(' ',1)
			name,genome,phage = seq_name.split(' ',1)[0].rsplit("_",2)
			XLS.write(name+'\t'+genome+'_'+phage+'\n')
			TMP.write('>'+name+"__"+genome+'_'+phage+' '+annotation+'\n'+s+'\n')
		TMP.close()
		shutil.move(TMP.name,outputname+".seq")
		TMP = open(outputname+'.tmp','w')
		for t,s in fasta_check(  open(outputname+'.pep'   )   ):
			seq_name,annotation = t[1:].strip().split(' ',1)
			name,genome,phage = seq_name.split(' ',1)[0].rsplit("_",2)
			TMP.write('>'+name+"__"+genome+'_'+phage+' '+annotation+'\n'+s+'\n')
		TMP.close()
		shutil.move(TMP.name,outputname+".pep")		
		
		
		con_data = Ddict()
		PHAGEXLS = open(outputname+"_PhageElement.xls",'w')
		PHAGEXLS.write("Name\tKind\tFunction\tRef_Source\tRef_Start\tRef_Stop\tRef_Frame\tSeq_Nucleotide\tSeq_Nucl_Length\n")
		for t,s in fasta_check( open(outputname+'.con','rU')):
			s = re.sub("\s+", "", s)
			[(start,end)] = re.findall( "\((\d+)\-(\d+) bp\)",t)
			length = int(end)-int(start)+1
			phage_name = t.split()[0].rsplit("_",1)[-1]
			phage_function = re.search("\S+\s+(.+)\s+\(",t).group(1)
			PHAGEXLS.write( prefix+'_'+phage_name+"\tPhageElement"+'\t'+phage_function+'\t'+prefix+'\t'+start+'\t'+end+'\t+\t'+s+'\t'+str(length)+'\n'  )
			
	else:
		README.write("""\n未发现任何前噬菌体序列。\n""")
			
	dir_list = glob.glob(outpath+'/*/')
	for e_dir in dir_list:
		shutil.rmtree(e_dir)
