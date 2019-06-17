#!/home/nfs/SOFTWARE/bin/nextflow

params.input = "./"
params.output = "./"
Result = params.output+"/Assembly/"



Channel.fromFilePairs(params.input+'/*.pair{1,2}').into {all_file }



process IDBA_ID_Denovo {
	executor 'pbs'
	cpus 32
	clusterOptions  "   -l nodes=1:ppn=16 -v PATH=$PATH "
	publishDir "$Result", mode: 'copy', overwrite: true
	input: 
		set val(sampleid),file(reads) from all_file
	output:
		file("*.contig.fa") into ass_result

		file("*.stats") into ass_stat
		
	script:
		"""	
		fq2fa --merge  ${reads[0]} ${reads[1]}  reads.fa
		idba_ud --read reads.fa  --pre_correction    -o  ${sampleid} --num_threads 32
		mv ${sampleid}/contig.fa ${sampleid}.contig.fa 
		N50 -i ${sampleid}.contig.fa  -o  ${sampleid}.stats
			 

		
		"""
		
}






