#!/usr/bin/nextflow
params.adapter = "/pub/SOFTWARE/Other/scythe/adap.fa"
params.quality = "20"
params.input = "./"
qc_path = params.input+"/Qc/"
stats_path = params.input+"/STATS/"


Channel.fromFilePairs(params.input+'/*.R{1,2}.*.gz').into {qc_plot_raw; all_file }



process RAW_QualityStats{
	
	maxForks 16
	input:
		set val(sampleid),file(reads) from qc_plot_raw
		
	output:
		file "*.dis" into raw_graph_stat
		file "*.stats" into  raw_stat
		
	"""
		zcat ${reads[0]} |  tee -a cache.R1.fastq | fastx_quality_stats -Q 33 -o ${sampleid}.1.dis
		zcat ${reads[1]} |  tee -a cache.R2.fastq | fastx_quality_stats -Q 33 -o ${sampleid}.2.dis


		Quality_Cat -1 cache.R1.fastq -2 cache.R2.fastq -n 33 -o  ${sampleid}.raw
		
		
	"""
	


}


process qc {
	maxForks 16
	publishDir "$qc_path", mode: 'copy', overwrite: true
	input: 
		set val(sampleid),file(reads) from all_file
	output:
		set val(sampleid),file("*pair{1,2}") into qc_result
		set val(sampleid),file("*pair{1,2}") into qc_result_ForStat
		
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

process QC_STAT_PROCESS{
	
	maxForks 16

	input:
	
		set val(sampleid),file(qc_pair) from qc_result_ForStat 
		
	output:

		file "*.stats" into qc_result_stat
		
		
		
	"""	

		Quality_Cat -1 ${qc_pair[0]} -2 ${qc_pair[1]} -n 33 -o  ${sampleid}.Filter
	"""



}





process QC_STAT_PROCESS{

	maxForks 16


	input:
	
		set val(sampleid),file(qc_pair) from qc_result 
		
	output:
		file "*.dis" into  qc_graph_stat

		
		
		
	"""	
		cat ${qc_pair[0]} | fastx_quality_stats -Q 33 -o ${sampleid}.Filter1.dis
		cat ${qc_pair[1]} | fastx_quality_stats -Q 33 -o ${sampleid}.Filter2.dis
		
	"""



}




		
filter_stat_graph = raw_graph_stat.mix(qc_graph_stat).flatten()



process QC_Plot {
	maxForks 16
	publishDir "$stats_path", mode: 'move', overwrite: true

	input:
		file filter_stats_fastx from filter_stat_graph
		
	output:
		file "*.pdf" into qc_graph

		
	"""
		fastx_quality_boxplot.R ${filter_stats_fastx} ${filter_stats_fastx}.qualstats.pdf&&fastx_nucleotide_distributionPer.R ${filter_stats_fastx} ${filter_stats_fastx}.nucdistr.pdf
		
		
	"""

		
}
total_stat = raw_stat.mix( qc_result_stat ).collect()
process Qc_Report{
	maxForks 16
	publishDir "$stats_path", mode: 'move', overwrite: true
	input:
		file all_statfile from total_stat
	
	output:
		file "*.tsv" into qc_report
		"""
			QC_ResultIntegrated.py  $all_statfile QC_Result.tsv
			
		"""
	}





