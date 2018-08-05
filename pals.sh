PALS_Align.py $1 ./pals_hit.gff3 >run.log 2>&1
	piler-64 -trs pals_hit.gff3 -out trs.gff  >run.log 2>&1
	mkdir Repeat 
	piler-64  -trs2fasta trs.gff -seq $1  -path Repeat >run.log 2>&1
	mkdir Mafft

	for i in Repeat/*.* ; do echo   "mafft  $i > Mafft/`basename $i.mafft`" >>mafft.sh; done
	cat mafft.sh| parallel -j 64 >run.log 2>&1
	mkdir Cons
	for i in Mafft/*.mafft; do piler-64 -cons $i -out Cons/`basename ${i%.*}.cons` -label `basename ${i%.*}`;done
    cat Cons/*.* >Total_Piler.fa
