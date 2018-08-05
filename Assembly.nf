#!/home/nfs/SOFTWARE/bin/nextflow

params.input = "./"
params.output = "./"
Result = params.output+"/Assembly/"



Channel.fromFilePairs(params.input+'/*.pair{1,2}').into {all_file }



process SoapDenovo {
	executor 'pbs'
	cpus 8
	clusterOptions  "   -l nodes=1:ppn=16 -v PATH=$PATH "
	publishDir "$Result", mode: 'copy', overwrite: true
	input: 
		set val(sampleid),file(reads) from all_file
	output:
		file("*.contig.fa") into ass_result

		file("*.stats") into ass_stat
		
	script:
		"""	
			 Multi_Soap.py -c 32 -t 32 -1  ${reads[0]}  -2   ${reads[1]} -n -i 400 -m 100 -r 80 -o Assembly 45,
			 mv   Assembly/K45/*.scafSeq ${sampleid}.contig.fa
			 mv Assembly/run*.log  ${sampleid}.stats
			 

		
		"""
		
}






