cnvnator  -root $1.root -tree $1.bam
 cnvnator  -root $1.root -his 100 -d $1
cnvnator -root $1cnvnator -root M2.root -stat 100    .root -stat 100
cnvnator -root $1.root   -partition 100  -ngc

cnvnator -root $1.root -call 100 -ngc > $1.CNV.out
 GFF_Region_Parse.py  $2 $1.Junction.Gt.vcf  $1.CNV.out
