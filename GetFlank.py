#!/usr/bin/python
from lpp import *
RAW = fasta_check(  open(sys.argv[1],'rU')  )
END = open("Flank.fa",'w')
all_length = {}
all_seq = {}
for t,s in RAW:
    seq = re.sub("\s+","",s)
    name = t.split()[0][1:]
    all_length[name]=len(seq)
    all_seq[ name ] = seq
    
for line in open(sys.argv[2],'rU'):
    qstart =0
    qend = 0
    line_l = line.strip().split("\t")
    length = int(line_l[4])
    start = int(line_l[2])-1
    end = int(line_l[3])-1
    gene  = line_l[-2]
    scaff = line_l[-1]
    if start >2000:
        qstart = start -2000
        start =2000
    if all_length[scaff] >(end+2000):
        qend = end+2000
    else:
        qend = all_length[scaff]
    end = start+length
    END.write(  '>'+gene+" from %s(%s bp) %s--%s   Gene location %s--%s\n"%(scaff,all_length[scaff],qstart,qend,start,end  )  )
    END.write(  all_seq[scaff][qstart:qend]+'\n')
