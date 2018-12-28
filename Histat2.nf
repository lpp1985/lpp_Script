#!/usr/bin/nextflow
params.input ="./"
Result_path = params.input+"/Result/"

db_name = file(params.index).name
db_path = file(params.index).parent
Channel.fromFilePairs(params.input+'/*R{1,2}*.gz').into { all_reads }


process mapping {
      executor 'pbs'
	maxForks 40
  cpus 32
  clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH"
  publishDir "$Result_path", mode: 'copy', overwrite: true
    input:
		set val(sampleid),file(reads) from all_reads

    output:
		file "*.bam" into STAR_Bam

    script:
		//
		// STAR Mapper
		//
		"""
		hisat2 -p 32  --dta -x $db_path/$db_name   -1  ${reads[0]} -2  ${reads[1]} -S align.sam
		samtools view -bS -F -@ 32 align.sam -o ${sampleid}.raw
		samtools sort  -@ 20 ${sampleid}.raw -o ${sampleid}.bam
		"""
   
}




process Combine{
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
	executor 'pbs'
	input:
		file all_bam from STAR_Bam.collect()
	output:
		file "Total.bam" into  bam_result
		file "Total.bam" into string_input
	script:
		"""
			samtools merge -@ 32  Total.bam $all_bam
			samtools sort  Total.bam  -o Total.sort.bam
			mv Total.sort.bam Total.bam
		"""


}

process StringTie{
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
	executor 'pbs'

	publishDir "$Result_path", mode: 'copy', overwrite: true
	input :
		file "total.bam" from string_input
		
	output:
		file "stringtie.gtf" into string_result

	script:
		"stringtie  total.bam   -o stringtie.gtf"


}

process BamHints{
	executor 'pbs'
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"

	publishDir "$Result_path", mode: 'copy', overwrite: true
	input :
		file "total.bam" from bam_result
		
	output:
		file "hints.gtf" into hints
		file "AlignStats.tsv" into stats
		
		
	script:
	
		"""
			bam2hints  --in=total.bam  --out=hints.gff
			bamtools stats -in  total.bam > AlignStats.tsv
		"""
		




}


