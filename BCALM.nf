#!/usr/bin/nextflow
params.input ="./"
params.transgene ="transgene.fasta"
Result_path = params.input+"/Result/"
genomeFile = file(params.transgene)
params.kmer ="41"
params.step = "6"
kmer= params.kmer
step = params.step
Channel.fromFilePairs(params.input+'/*_{1,2}.fq.gz').into { all_reads }


process BCALM {
	publishDir "$Result_path", mode: 'copy', overwrite: true
    	executor 'pbs'
	cpus 32
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"

    input:
		
		set val(sampleid),file(reads) from all_reads


    output:
		file "list_reads*.fa" into Graph

    script:

		"""
		ls  ${reads[0]}  ${reads[1]}  >list_reads
		bcalm -kmer-size $kmer -nb-cores 64 -in list_reads -abundance-min 11 -max-memory 100000
		



"""
   
}

process Creep{
	publishDir "$Result_path", mode: 'move', overwrite: true
	executor 'pbs'
 	cpus 1
 	clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH"

        input:
                file all_graph from Graph
				file genomeFile
				val kmer from params.kmer
				val step from params.step = "6"
        output:
                file "Graph.fa" into result
        script:
                """
					ContainRef   -i $all_graph  -l $genomeFile -o test.list 
					Kmer_Graph -i $all_graph -k $kmer -l test.list -o Graph.fa -s $step
                """


}






