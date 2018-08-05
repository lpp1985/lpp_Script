#!/home/nfs/SOFTWARE/bin/nextflow
params.query = "/sample.fa"

params.out = "./genblasta.out"
params.db ="Assem"
db_name = file(params.db).name
db_path = file(params.db).parent
 
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
Channel
    .fromPath(params.query)
    .splitFasta(by: 500)
    .set { fasta }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */
process GenBlastA {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
		file 'query.fa' from fasta
		
 
    output:
		  file 'top_hits' into top_hits
 
    """
		genblasta -P blast -pg tblastn -q query.fa  -t Assem -o top_hits
    
    """
}
top_hits.collectFile(name: params.out)

