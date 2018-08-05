#!/home/nfs/SOFTWARE/bin/nextflow
params.Input = "./*.fasta"

params.out = "./Cazy"
params.eval = "1e-5"


Channel.fromPath(params.Input).into {all_protein }

process Nr_Annotation {
    executor 'pbs'
	publishDir "${params.out}", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
    input:
		file protein from all_protein
	output:
		file "*.annot" into Result

	script:
		out_name = protein.baseName
    """
	  hmmscan --cpu 16  --domtblout ${out_name}.dbCAN -E 1e-5 /home/nfs/Database/Cazy/dbCAN-fam-HMMs.txt $protein  >/dev/null
	  hmmscan-parser.sh ${out_name}.dbCAN >${out_name}.annot
    
    """
}
