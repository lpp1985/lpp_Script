.  /home/nfs/SOFTWARE/Pacbio/smrtanalysis/install/smrtanalysis_2.3.0.140936/etc/setup.sh 
for i in *.bax.h5; do bash5tools.py    --outType   fasta --minReadScore 0.80 --outFilePrefix ${i%.*.*} $i&done
