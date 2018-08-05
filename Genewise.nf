#!/home/nfs/SOFTWARE/bin/nextflow
params.command = "$HOME/sample.fa"
params.out = "./result.txt"


 
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
Channel
    .fromPath(  params.command   )splitText(by: 10)
    .set { command }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */
process genewise {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH"
    input:
		file "run.sh" from command
		
 
    output:
		 stdout channel
 
    """
	bash run.sh
    
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
channel.collectFile(name: params.out)

