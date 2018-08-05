#!/home/nfs/SOFTWARE/bin/nextflow
params.Input = "./*.bam"
params.Bed = "unigene.bed"
params.out = "./Out"



Channel.fromPath(params.Input).into {bamall }
bedall=file(params.Bed)


process RSEQC{
	executor 'pbs'
	publishDir "${params.out}/Rseqc", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
	input:
		file bam from bamall
		file bed  from bedall
	output:
		file "*.pdf" 
		file "*.stats"
		file "*.geneBodyCoverage.r" into genebodyr
	script:
		out_name = bam.baseName
	 """
	  samtools index ${bam}
	  inner_distance.py  -i $bam -r $bed -o ${out_name}.Inner 
	  RPKM_saturation.py  -i $bam -r $bed  -o ${out_name}.saturation 
	  geneBody_coverage.py  -i $bam  -r $bed -o ${out_name}
	  bamtools stats -in $bam > ${out_name}.stats
	 """

}
all_r_file = genebodyr.collect()
process GeneBody{
	executor 'pbs'
	publishDir "${params.out}/Rseqc", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH"
	input:
		file r_file from all_r_file
	output:
		file "*.pdf" 

	script:
	 """
	  geneBody_coveragePlot.py $r_file Total
	  
	 """

}


