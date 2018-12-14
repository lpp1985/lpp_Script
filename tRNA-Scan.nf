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
process tRNA {
    executor 'local'
    cpus 8
    clusterOptions  " -d $PWD  -l nodes=1:ppn=2 -v PATH=$PATH"
    input:
    file 'query.fa' from fasta
    file db_path
 
    output:
    file 'tRNA.out' into top_hits

 
    """
     tRNAscan-SE -G -o tRNA.out -f tRNA.structure -m stats -q  query.fa        
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
process Parse{
	publishDir "$output_path", mode: 'copy', overwrite: true
	executor 'pbs'
	cpus 2
	clusterOptions  " -d $PWD  -l nodes=1:ppn=2 -v PATH=$PATH"
	input:
		file "tRNAScan.out" from Result
	output:
		file "*.gff" into taxon

	"tRNAscan-SE_Gff.py  tRNAScan.out  tRNA.gff"
}
