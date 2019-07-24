#!/usr/bin/nextflow
params.input ="./"
params.genome ="scaff.fa"
params.index = "REF*"
params.output= "Result"
Result_path = params.output
genomeFile = file(params.genome)
db_name = file(params.index).name
db_path = file(params.index).parent
Channel.fromFilePairs(params.input+'/*_R{1,2}.gz').into { all_reads; raw_reads }


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
		bwa mem  -M  -t 64  $db_path/$db_name   ${reads[0]}  ${reads[1]} 1>${sampleid}.raw  2>/dev/null
		
		samtools view -bS -@ 20  ${sampleid}.sam   -o  ${sampleid}.raw   2>/dev/null
		samtools sort ${sampleid}.raw  -o ${sampleid}.bam
		HTseq_Count.py -s  ${sampleid}.bam -o ${sampleid}.count
		PickBestRef.py  ${genomeFile}  ${sampleid}.count



"""
   
}

process Polish_mapping {
    	executor 'pbs'
	cpus 16
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
	publishDir "$Result_path", mode: 'copy', overwrite: true

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
		Polishing.sh ${ref} ${sampleid}

"""
   
}





