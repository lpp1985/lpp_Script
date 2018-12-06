#!/home/nfs/SOFTWARE/bin/nextflow 
params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/pdb/pdb"
params.out = "./result.txt"
params.chunkSize = 5000
db_name = file(params.db).name
db_path = file(params.db).parent
output_path = file(params.out).parent+"/"
Result = Channel.create() 
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
Channel
    .fromPath(params.query)
    .splitFasta(by: params.chunkSize)
    .set { fasta }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */
process blast {
    executor 'local'
    cpus 8
    clusterOptions  " -d $PWD  -l nodes=1:ppn=2 -v PATH=$PATH"
    input:
    file 'query.fa' from fasta
    file db_path
 
    output:
    file 'top_hits' into top_hits
 
    """
    /usr/bin/blastn -task 	blastn    -query query.fa -db $db_path/$db_name   -num_threads  10 -outfmt 6 -max_target_seqs 1   -out  top_hits -evalue 1e-35
    
    """
}
 
 
/*
 * Each time a file emitted by the 'top_hits' channel an extract job is executed
 * producing a file containing the matching sequences
 */

/*
 * Collects all the sequences files into a single file
 * and prints the resulting file content when complete
 */
top_hits.collectFile(name: params.out).into(Result)
process Taxon{
	publishDir "$output_path", mode: 'copy', overwrite: true
	executor 'pbs'
	cpus 2
	clusterOptions  " -d $PWD  -l nodes=1:ppn=2 -v PATH=$PATH"
	input:
		file "Blast.out" from Result
	output:
		file "*.tax" into taxon

	"/home/nfs/SOFTWARE/Other/assemblage/blast_taxonomy_report.pl -b Blast.out    -nodes /home/nfs/Database/Taxon/nodes.dmp    -names  /home/nfs/Database/Taxon/names.dmp   -gi_taxid_file /home/nfs/Database/Taxon/gi_taxid_nucl.dmp    -t genus=1 -t order=1 -t family=1 -t superfamily=1 -t kingdom=1  >megablast.nt.tax "


}
