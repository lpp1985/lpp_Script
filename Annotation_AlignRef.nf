#!/usr/bin/nextflow
params.query = "./"
params.ref="ref.fasta"
params.output="./Annotation_And_Ref/"
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
outputpath = params.output
ref_data = file( params.ref)

Channel
    .fromPath(params.query)
    .into { Prokka_input }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */

 process PROKKA {
	maxForks 40
	publishDir "${outputpath}/", mode: 'copy', overwrite: true
	executor 'pbs'
        cpus 1
        clusterOptions  "   -l nodes=1:ppn=32 -v PATH=$PATH;PERLLIB=$PERLLIB "
    input:

		file Scaff from Prokka_input
		file ref_data
 
    output:
		file '*.tar.gz' into Annotationa
    script:
	  res_path=Scaff.getBaseName()
 
    """	mkdir ${res_path}
	
	
	prokka   $Scaff --force  --prefix $res_path --outdir $res_path --evalue 1e-5  -cpus 64 --quiet  
    /home/nfs/SOFTWARE/Other/MicroPipe/Scripts/create_database.py    -n ${res_path}/*.ffn -p ${res_path}/*.faa -f ${res_path}/*.fna -g ${res_path}/*.gff -d ${res_path}/Gene.xls

	mkdir ${res_path}/REF

	promer --maxmatch $ref_data $Scaff
	show-coords -oTl out.delta >${res_path}/REF/${res_path}.align.xls
	show-coords -oTH out.delta | grep  '\\['  |cut -f13 >Need.tsv
	promer --maxmatch $ref_data ${res_path}/*.ffn 
        show-coords -oTl out.delta >${res_path}/REF/${res_path}.cds.xls
	
	
	GetSeq -f $Scaff -d Need.tsv -o ${res_path}/REF/${res_path}.Ref.fa -n 1

	GetLine.py ${res_path}/Gene.xls Need.tsv ${res_path}/REF/Annotation.xls
	tar -zcf ${res_path}.tar.gz ${res_path}/
    """
}

 
 
 
