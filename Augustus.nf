#!/usr/bin/nextflow
params.input ="./"

Result_path = params.input+"/GenePrediction/"


all_genome= file(params.input+'/*.fna')


process Augustus {
		executor 'local'
		maxForks 4
		cpus 32
 
		publishDir "$Result_path", mode: 'copy', overwrite: true
    input:
		
		file ref from all_genome

    output:
		file "*.gffo" into GFF
		file "*.cds" into CDS
		file "*.pep" into PEP

    script:
		out_name = ref.baseName
		
		"""
		mkdir  Ref
		splitMfasta.pl $ref  --outputpath=Ref --minsize=1
		cd Ref
		for i in *.fa; do echo "/home/nfs/SOFTWARE/Other/augustus/augustus/bin/augustus  --AUGUSTUS_CONFIG_PATH=/home/nfs/SOFTWARE/Other/augustus/augustus/config/ --species=TCK  \$i 1>\$i.gff">>run.sh ; done
		 cat run.sh  | parallel -j 10
		 cat *.gff >augustus.tmp.gff3
		join_aug_pred.pl  < augustus.tmp.gff3 >augustus.gff
		gtf2gff.pl  <augustus.gff  --gff3 --printExon  --out ${out_name}.gffo
		mv ${out_name}.gffo ../
		cd ../
		
		gffread  -g $ref  -x ${out_name}.cds ${out_name}.gffo
		gffread  -g $ref  -x ${out_name}.pep ${out_name}.gffo
		"""
   
}






