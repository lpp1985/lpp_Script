#!/usr/bin/nextflow
params.input ="./"
params.genome ="scaff.fa"
Result_path = params.input+"/Result/"
genomeFile = file(params.genome)

Channel.fromFilePairs(params.input+'/*_1{1,2}*.gz').into { all_reads }
process index {
    executor 'pbs'
    cpus 32 
    clusterOptions  "   -d $PWD "
    input:
		file genomeFile

    
    output:
		file "STARgenome" into STARgenomeIndex

    script:
    //
    // STAR Generate Index and create genome length file
    //
    """
        mkdir STARgenome
        STAR --runThreadN 21 \
             --runMode genomeGenerate \
             --genomeDir STARgenome \
             --genomeFastaFiles ${genomeFile} \
             --outFileNamePrefix STARgenome

    """
}


process mapping {
      executor 'pbs'
	maxForks 4
  cpus 32
  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -V"
  publishDir "$Result_path", mode: 'copy', overwrite: true
    input:
		file STARgenome from STARgenomeIndex.first()
		set val(sampleid),file(reads) from all_reads

    output:
		file "*.bam" into STAR_Bam

    script:
		//
		// STAR Mapper
		//
		"""
		STAR  --runThreadN 20 --readFilesCommand   zcat  --readFilesIn ${reads[0]} ${reads[1]}  --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 10000 --alignMatesGapMax 10000 --chimSegmentMin 20 --twopassMode Basic --outFileNamePrefix ./Star_OUT \
		--genomeDir $STARgenome --outSAMstrandField intronMotif RemoveNoncanonical --alignTranscriptsPerReadNmax 10000  --limitBAMsortRAM 8921722571 --outSAMunmapped Within
		mv *.bam  ${sampleid}.bam
		"""
   
}




process Combine{
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -V"
	executor 'pbs'
	input:
		file all_bam from STAR_Bam.collect()
	output:
		file "Total.bam" into  bam_result
		file "Total.bam" into string_input
	script:
		"""
			samtools merge -@ 32  Total.bam $all_bam
			/home/nfs/SOFTWARE/bin/samtools sort -@ 32 -m 4G Total.bam  -o Total.sort.bam
			mv Total.sort.bam Total.bam
		"""


}

process StringTie{
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -V"
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
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -V"

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


