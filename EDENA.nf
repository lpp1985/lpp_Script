#!/usr/bin/nextflow
params.input ="./"

Result_path = params.input+"/Result/"


Channel.fromFilePairs(params.input+'/*_R{1,2}.fq').into { all_reads }



process EdenaAssem {
	publishDir "$Result_path", mode: 'copy', overwrite: true
    	executor 'pbs'
	cpus 32
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"

    input:

		set val(sampleid),file(reads) from all_reads

    output:
		file "*.fasta" into total_bam

    script:

		"""
		edena  -paired  ${reads[0]}  ${reads[1]} -p ${sampleid} -nThreads 16
		edena -cc yes  -e *.ovl -p ${sampleid}
		mv *.fasta ${sampleid}.fasta
		


"""
   
}

