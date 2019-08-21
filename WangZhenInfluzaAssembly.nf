#!/usr/bin/nextflow
params.input ="./"
params.genome ="scaff.fa"
params.index = "REF*"
params.output= "./Result"
Result_path = params.output
genomeFile = file(params.genome)
db_name = file(params.index).name
db_path = file(params.index).parent
Channel.fromFilePairs(params.input+'/*_{1,2}.fq.gz').into { all_reads; raw_reads }


process raw_mapping {
    	executor 'pbs'
	cpus 32
	clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH"

    input:
		file genomeFile
		set val(sampleid),file(reads) from all_reads


    output:
		file "*.bam" into bam

    script:

		"""
		bwa mem  -M  -t 32  $db_path/$db_name   ${reads[0]}  ${reads[1]} 1>${sampleid}.sam  2>/dev/null
		
		samtools view -bS -@ 20  ${sampleid}.sam   -o  ${sampleid}.raw   2>/dev/null
		samtools sort ${sampleid}.raw  -o ${sampleid}.bam



"""
   
}
process Pick_Best{
	executor 'pbs'
	cpus 1
	clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH"

    input:
		file genomeFile	
		file raw_bam from bam
	output:
		file "Best1.fa" into best_ref


	script:
		sample = raw_bam.baseName
		"""
		
		HTseq_Count.py -s  ${raw_bam} -o ${sample}.count
		PickBestRef.py  ${genomeFile}  ${sample}.count



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
		sed -ir "s/_pilon//g" *.fasta

"""
   
}





