#!/usr/bin/nextflow
params.input ="./"
params.genome ="scaff.fa"
Result_path = params.input+"/Result/"
genomeFile = file(params.genome)

Channel.fromFilePairs(params.input+'/*_{1,2}.fq.gz').into { all_reads }
process index {
    executor 'local'
    scratch false
    cpus 40 
    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V "
    input:
		file genomeFile

    
    output:
		file "*.bt2*" into STARgenomeIndex

    script:
    //
    // Bowtie2 Generate Index and create genome length file
    //
    """
	  bowtie2-build --threads  40 ${genomeFile} REF
     

    """
}


process mapping {
	publishDir "$Result_path", mode: 'copy', overwrite: true
    	executor 'local'
	maxForks 4
	cpus 4
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
 
    input:
		file genomeFile
		file STARgenome from STARgenomeIndex
		set val(sampleid),file(reads) from all_reads

    output:
		file "*.bam" into total_bam
		file "*.stats" into stats

    script:

		"""
		bowtie2 -X 1000 --end-to-end --very-sensitive  --mm --reorder -t -p 40 -x  REF -1 ${reads[0]} -2  ${reads[1]} |samtools view -bS -F12 -@40 -T ${genomeFile} - | samtools sort -@10 -o ${sampleid}.bam 
		bamtools stats -in ${sampleid}.bam > ${sampleid}.stats
		



"""
   
}

