#!/usr/bin/nextflow
params.input ="./"
params.genome ="scaff.fa"
Result_path = params.input+"/Result/"
genomeFile = file(params.genome)

Channel.fromFilePairs(params.input+'/*_{1,2}.fq.gz').into { all_reads, raw_reads }


process raw_mapping {

    executor 'pbs'
	cpus 32
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"

    input:
		file genomeFile
		
		set val(sampleid),file(reads) from all_reads
		file genomeFile

    output:
		file "Best1.fa" into best_ref

    script:

		"""
		bwa index  ${genomeFile} -p REF
		bwa mem  -M  -t 64  REF   ${reads[0]}  ${reads[1]}  1> ${sampleid}.sam  2>/dev/null
		
		samtools view -bS -@ 20  ${sampleid}.sam   -o  ${sampleid}.raw   2>/dev/null
		samtools sort ${sampleid}.raw  -o ${sampleid}.bam
		HTseq_Count.py -i  ${sampleid}.bam -o ${sampleid}.count
		PickBestRef.py  ${genomeFile}  ${sampleid}.count



"""
   
}
process index_pol {
    executor 'pbs'
    scratch true
    cpus 1 
    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V "
    input:
		file "Best1.fa" from  best_ref

    
    output:
		file "REF*" into STARgenomeIndex

    script:
    //
    // STAR Generate Index and create genome length file
    //
    """
        bwa index  Best1.fa -p REF

    """
}

process Polish_mapping {
    	executor 'pbs'
	cpus 16
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"

    input:
		file ref from  best_ref
		
		set val(sampleid),file(reads) from raw_reads

    output:
		file "*.fasta" into pol_bam

    script:

		"""
		bwa index  $ref -p REF
		bwa mem  -M  -t 16  REF   ${reads[0]}  ${reads[1]}  1> ${sampleid}.sam  2>/dev/null
		
		samtools view -bS -@ 20  ${sampleid}.sam   -o  ${sampleid}.raw   2>/dev/null

		samtools sort ${sampleid}.raw  -o ${sampleid}.bam
		samtools index  ${sampleid}.bam
		Polishing.sh ${genomeFile} ${sampleid}

"""
   
}





