#!/home/nfs/SOFTWARE/bin/nextflow
params.command = "$HOME/sample.fa"
params.out = "./result.txt"
/*output_path= file(params.command).parent+"/"
result = output_path+'/'+params.out 
*/
 
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
Channel
    .fromPath(  params.command   )splitText(by: 100)
    .set { command }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */
process PALS {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=2 -v PATH=$PATH"
    input:
		file "run.sh" from command
		
 
    output:
		 file "result.tsv" into top_hits
 
    """
	bash run.sh
	cat *.xml > result.tsv
    
    """
}
 
top_hits.collectFile(name: params.out) 

