#!/home/nfs/SOFTWARE/bin/nextflow
params.Protein = "./*.fasta"
params.Nucl = "./*.fna"
params.out = "./Out"
params.eval = "1e-5"


Channel.fromPath(params.Protein).into {protein_Cazy; protein_TMHMM; protein_VFDB;  protein_CARDB; protein_PHI }
Channel.fromPath(params.Nucl).into {nucl_Resfinder}
process Cazy {
    executor 'pbs'
	publishDir "${params.out}/Cazy", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
    input:
		file protein from protein_Cazy
	output:
		file "*.annot" into CazyResult

	script:
		out_name = protein.baseName
    """
	  hmmscan --cpu 16  --domtblout ${out_name}.dbCAN -E 1e-5 /home/nfs/Database/Cazy/dbCAN-fam-HMMs.txt $protein  >/dev/null
	  hmmscan-parser.sh ${out_name}.dbCAN >${out_name}.annot
    
    """
}
process TMHMM {
    executor 'pbs'
	publishDir "${params.out}/TMHMM", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
    input:
		file protein from protein_TMHMM
	output:
		file "*.TMHMM.tsv" into TMHMMResult

	script:
		out_name = protein.baseName
    """
	 THMMPrediction.py  -i $protein -o ${out_name}.TMHMM.tsv 
    
    """
}
process VFDB{
	executor 'pbs'
	publishDir "${params.out}/VFDB", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
	input:
		file protein from protein_VFDB
	output:
		file "*.VFDB.Result" into VFDBResult
	script:
		out_name = protein.baseName
	 """
	 echo -e "Name\tVFDBID\tVFDB_Annotation\tVFDB_Identity\tVFDB_Evalue\n" >${out_name}.VFDB.Result
	 blastp  -query $protein -db /home/nfs/Database/VFDB/VFDB  -evalue 1e-35 -outfmt  "6 qseqid sacc  stitle pident  evalue" -num_threads 12  -num_alignments  1  >>${out_name}.VFDB.Result
	 """

}

process CARDB {
    executor 'pbs'
	publishDir "${params.out}/CARDB", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
    input:
		file protein from protein_CARDB
	output:
		file "*.CARDB.txt" into CARDBResult

	script:
		out_name = protein.baseName
    """
	 rgi.py  -t protein -i $protein  -n 64  -o  ${out_name}.CARDB  -a blast
    
    """
}


process PHI{
	executor 'pbs'
	publishDir "${params.out}/PHI", mode: 'copy', overwrite: true
    clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
	input:
		file protein from protein_PHI
	output:
		file "*.PHI.tsv" into PHIBResult
	script:
		out_name = protein.baseName
	 """
	  PHI.py  -i $protein   -o ${out_name}.PHI.tsv -x /home/nfs/Database/PHI/phi_accessions.txt 
	  
	 """

}
process ResFinder{
	executor 'pbs'
	publishDir "${params.out}/Resfinder", mode: 'copy', overwrite: true
	clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH"
	input :
		file nucl from nucl_Resfinder
	output:
		file "*.tar.gz" into ResFinerResult

	script:
		out_name = nucl.baseName
		""" resfinder.py -i $nucl  -o ./ -p /home/nfs/Database/resfinder_db/
		    mkdir $out_name
		    mv *_*  ${out_name}/
		    tar -zcf ${out_name}.tar.gz  ${out_name}/
		"""
}
