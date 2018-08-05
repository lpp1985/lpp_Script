#!/usr/bin/nextflow
params.query = "./"
params.ref="ref.fasta"
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
outputpath = "./"
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
	publishDir "$outputpath/Annotation/", mode: 'copy', overwrite: true
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
 
    """
	 Contig_Combine.py  $Scaff  Scaff.fa
	 prokka   Scaff.fa  --force  --prefix $res_path --outdir $res_path --evalue 1e-5  --genus test --strain $res_path --cpus 64 --compliant --quiet  --locustag $res_path
	/home/nfs/SOFTWARE/Other/MicroPipe/Scripts/create_database.py    -n ${res_path}/*.ffn -p ${res_path}/*.faa -f ${res_path}/*.fna -g ${res_path}/*.gff -d ${res_path}/Gene.xls
	 nucmer --maxmatch $ref_data ${res_path}/*.ffn
 	show-coords -oTl out.delta >${res_path}/Ref_gene.tsv
	nucmer --maxmatch $ref_data ${res_path}/*.fna
	show-coords -oTl out.delta >${res_path}/Ref_Genome.tsv
	AnnotationPipe.py  -p ${res_path}/*.faa -n  ${res_path}/*.ffn   -o ${res_path}/Annotation -e 1e-5 
	XLS_Combine.py -i ${res_path}/Annotation/Detail/ -o ${res_path}/Annotation/Table -g ${res_path}/Gene.xls
	 tar -zcf ${res_path}.tar.gz $res_path/
	 
    """
}

 
 
 
