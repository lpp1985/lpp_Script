#!/usr/bin/nextflow
params.query = "$HOME/sample.fa"
params.out = "./result.txt"

Channel
    .fromPath(params.query)
    .into { RAW }
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */
process KaKs {
    executor 'pbs'
    cpus 1
    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH"
    input:
    file axt from RAW
 
    output:
        file '*.result' into top_hits
    script:
    """
	KaKs_Calculator  -i $axt -o test.axt -m NG 1>/dev/null ; cut -f 1,2,3,4,5 test.axt  >out.result
    
    """
}
top_hits.collectFile(name: params.out)
