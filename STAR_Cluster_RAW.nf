#!/usr/bin/nextflow
params.input ="./"
params.db ="scaff.fa"
Result_path = params.input+"/Result/"
dbFile = params.db

Channel.fromFilePairs(params.input+'/*_R{1,2}.fq.gz').into { all_reads }


process mapping {
	executor 'pbs'
	cpus 32
 	clusterOptions  "   -l nodes=1:ppn=32 -v PATH=$PATH "


    input:
		val dbFile 
		set val(sampleid),file(reads) from all_reads

    output:
		file "*.bam" into STAR_Bam

    script:
		//
		// STAR Mapper
		//
		"""
		STAR --runThreadN 41 --readFilesCommand   zcat  --readFilesIn ${reads[0]} ${reads[1]}  --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 --twopassMode Basic --outFileNamePrefix ./Star_OUT \
		--genomeDir $dbFile --outSAMstrandField intronMotif RemoveNoncanonical --alignTranscriptsPerReadNmax 100000  --limitBAMsortRAM 8921722571
		mv *.bam  ${sampleid}.bam
		"""
   
}



total_bam = STAR_Bam.collect()

process Combine{
 	executor 'pbs'
 	cpus 32
 	clusterOptions  "   -l nodes=1:ppn=32 -v PATH=$PATH "



	input:
		file all_bam from total_bam
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
	executor 'pbs'
	cpus 32
	clusterOptions  "   -l nodes=1:ppn=32 -v PATH=$PATH "

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
 	cpus 2
 	clusterOptions  "   -l nodes=1:ppn=2 -v PATH=$PATH "
	publishDir "$Result_path", mode: 'copy', overwrite: true
	input :
		file "total.bam" from bam_result
		
	output:
		file "hints.gtf" into hints
		file "Align.stats" into stats
		
		
	script:
	
		"""
			bam2hints  --in=total.bam  --out=hints.gff
			bamtools stats -in  total.bam > AlignStats.tsv
		"""
		




}



