#!/home/nfs/SOFTWARE/bin/nextflow

params.input = "./"
params.output = "./"
Result = params.output+"/MLST/"



 Channel.fromPath( '*.contig' ).into {all_file }



process MLST {
	executor 'local'
	cpus 8
	publishDir "$Result", mode: 'copy', overwrite: true
	maxForks 4
	input: 
		file each_Ctg from all_file
	output:
		file("*.html") into MLST_Result


		
	script:
		name = file(each_Ctg).name
		"""	
			Contig_Combine.py  $each_Ctg  Scaff.fa
			MLST.py -i  Scaff.fa -o ${name}
			rename "s/contig.//g" *.html
			
			 
			 
			 

		
		"""
		
}






