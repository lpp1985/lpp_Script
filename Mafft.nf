#!/home/nfs/SOFTWARE/bin/nextflow
params.Input = "./F*"



Channel.fromPath(params.Input,type: 'dir' ).into {all_path }

process Mafft {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH"
    input:
		file path from all_path


	script:
		
    """
		mafft ${path}/seq.fa >${path}/mafft
		trimal  -strict  -in  ${path}/mafft   -out ${path}/out   
    
    """
}
