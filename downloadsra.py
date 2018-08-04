#!/usr/bin/env python
from lpp import *
all_need = sys.argv[1]

commandline = "ascp -QT -i /home/nfs/SOFTWARE/bin/aspera/etc/asperaweb_id_dsa.openssh  anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra  . &&  for i in *.sra; do fastq-dump --split-3 `pwd`/$i& done\n"%( all_need[:3],all_need[:6],all_need,all_need  )
COM = open("download.sh",'a')
COM.write(commandline)
