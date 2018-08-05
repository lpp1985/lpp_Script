#!/home/nfs/SOFTWARE/bin/nextflow
params.query = "$HOME/sample.fa"
params.db = "$HOME/tools/blast-db/pdb/pdb"
params.out = "./result.txt"
params.chunkSize = 500
db_name = file(params.db).name
db_path = file(params.db).parent
 
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
    executor 'pbs'
    cpus 2 
    clusterOptions  " -d $PWD   -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
    file 'query.fa' from fasta
 
    output:
    file 'top_hits' into top_hits
 
    """
	 genblasta -P blast -pg tblastn -q query.fa -p T  -t $db_path/$db_name -o top_hits
    
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
top_hits.collectFile(name: params.out)

