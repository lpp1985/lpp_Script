#!/home/nfs/SOFTWARE/bin/nextflow
params.adapter = "/home/nfs/SOFTWARE/Other/scythe/illumina_adapters.fa"
params.quality = "20"
params.input = "./"
qc_path = params.input+"/Qc/"
stats_path = params.input+"/STATS/"


Channel.fromFilePairs(params.input+'/*_{1,2}.*.gz').into { all_file }





process qc {
	maxForks 16
	
	input: 
		set val(sampleid),file(reads) from all_file
	output:
		set val(sampleid),file("*pair{1,2}") into qc_result_gz
	script:
		"""
		sickle pe -t sanger -n -l 100 -q ${params.quality} \
			-f <(scythe -a ${params.adapter} <(zcat ${reads[0]}) -q sanger 2> /dev/null )\
			-r <(scythe -a ${params.adapter} <(zcat ${reads[1]}) -q sanger 2> /dev/null )\
			-o  ${sampleid}.pair1\
			-p  ${sampleid}.pair2\
			-s  ${sampleid}.fail

		"""
		
}

process GZ{
	
	maxForks 20
	publishDir "$qc_path", mode: 'copy', overwrite: true

	input:
	
		set val(sampleid),file(qc_pair) from qc_result_gz 
		
	output:

		file ("*R?.fq.gz") into gzip_result
		
		
		
	"""	
		gzip -c ${qc_pair[0]}  >${sampleid}_R1.fq.gz

		gzip -c ${qc_pair[1]} > ${sampleid}_R2.fq.gz

		
		
	"""



}



