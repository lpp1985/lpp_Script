from lpp import *
RAW = fasta_check( open(sys.argv[1])  )
all_seqhash = {}
for t,s in RAW:
	all_seqhash[t[1:].split()[0]]=s

RAW = open(sys.argv[2])
for line in RAW:
	line_l = line.split("\t")
	name = line_l[0]
	if not os.path.exists(name):
		os.mkdir(name)
	path = name+'/'
	END = open(path+"/seq.fa",'w')
	all_seq = sorted(line_l[1].split(),key=lambda x: x.split('|')[0])
	for seq in all_seq:
		END.write(">"+seq.split("|")[0]+'\n'+ all_seqhash[seq]  )

