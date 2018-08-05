#!/usr/bin/nextflow
params.query = "./"
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
outputpath = "./"

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
        clusterOptions  "   -l nodes=1:ppn=32 -V "
    input:

		file Scaff from Prokka_input
 
    output:
		file '*.tar.gz' into Annotationa
    script:
	  res_path=Scaff.getBaseName()
 
    """
	 Contig_Combine.py  $Scaff  Scaff.fa
	export PERL5LIB=$PERL5LIB:/home/nfs/SOFTWARE/Other/prokka-1.11/bin
	 prokka   Scaff.fa  --force  --prefix $res_path --outdir $res_path --evalue 1e-5  --genus test --strain $res_path --cpus 64 --compliant --quiet  --locustag $res_path
	/home/nfs/SOFTWARE/Other/MicroPipe/Scripts/create_database.py    -n ${res_path}/*.ffn -p ${res_path}/*.faa -f ${res_path}/*.fna -g ${res_path}/*.gff -d ${res_path}/Gene.xls

 	
	AnnotationPipe_NOSPHIX.py  -p ${res_path}/*.faa -n  ${res_path}/*.ffn   -o ${res_path}/Annotation -e 1e-5 
	XLS_Combine.py -i ${res_path}/Annotation/Detail/ -o ${res_path}/Annotation/Table -g ${res_path}/Gene.xls
	PGAP_Prepare.py   ${res_path}
	mv *.pep ${res_path}
	mv *.function ${res_path}
	mv *.nuc ${res_path}
	tar -zcf ${res_path}.tar.gz $res_path/
	 
    """
}

 
 
 
