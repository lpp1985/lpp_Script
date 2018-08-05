#!/usr/bin/nextflow
params.query = "$HOME/Scaffold.fa"
params.repmask_lib = "/home/nfs/SOFTWARE/Other/RepeatMasker/Libraries/RepeatMasker.lib"
 
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
scaff = file("$params.query")
outputpath = scaff.getParent()

 Channel
    .fromPath(params.query)
    .into { scaffold_pro; scaffold_build ;scaffold_den;scaffold_raw; scaffold_ltr; scaffold_piler }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */

 process PILER_DP {
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V"

    input:

		file 'Scaffolds.fa' from scaffold_piler
 
    output:
		file 'Total_Piler.fa' into piler_lib
 
    """
	pals.sh Scaffolds.fa
   
    """
}

 
 
 
 
process LTR_FINDER {
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH;PERLLIB=$PERLLIB"


    input:

		file 'scaff.fa' from scaffold_ltr
 
    output:
		file '*.gff3' into ltr_gff
 
    """
    ltr_finder -w 2 scaff.fa  > ltr_raw
	LTR_PARSE.py -i  ltr_raw  -o  ltr.gff3
    
    """
}



process RepeatPro {
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -V"


    input:
		file 'scaff.fa' from scaffold_pro

 
    output:
		file '*.gff3' into repeatmaker_pro
 
    """
    RepeatProteinMask      -engine  ncbi scaff.fa
	RepeatProteinMaskResultTrans.py  *.annot Protein.gff3
    
    """
}

process RepeatDenovo_BuildDatabase {
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V"


    input:
		file 'scaff.fa' from scaffold_build

 
    output:
		file 'Database*' into database
 
    """
     BuildDatabase -name Database  scaff.fa
    
    """
}

process RepeatDenovo_Modeler{
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V"
	input :
		file db from database
	output :
		file "RM_*/consensi.fa.classified" into denovobase
	"""
	 RepeatModeler -database Database
	"""

}

process RepeatDb_IntegrateAnnotation{
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V"
	input :
		file "denovo.lib" from denovobase
		file  "piler.lib " from piler_lib
	output :
		file "RepBase.fa" into replib
	"""
	cat_all_fasta.py  -i ./ -a lib  -o Rep -n Rep
	cdhit-est -i Rep.fasta -M 0 -c 0.8 -o Clustered

	RepAnnotation.py  -i  Clustered -d  $params.repmask_lib  -o   RepBase.fa  						 
	"""

}


 
process RepeatDenovo_Masker{
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V"
	input:
		file "lib" from replib
		file 'scaff.fa' from scaffold_den
	output:
		file "Denovo/*.gff" into repeatmaker_deno
		file "Repeat.xls" into repeatmaker_table
	"""
		mkdir Denovo ; RepeatMasker -e ncbi  -pa 8 -lib lib -dir Denovo -gff scaff.fa
		RepMaskerr_Stats.py  -s scaff.fa  -i Denovo/*.out  -o Repeat.xls
	
	"""
		}



		
process Integrate{
	executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -V"

	publishDir "$outputpath/RepeatMasked/", mode: 'copy', overwrite: true
	input:
		file "Denovo.gff3" from repeatmaker_deno
		file "Protein.gff3" from repeatmaker_pro
		file "Scaffolds.fa" from scaffold_raw
		file "RepeatMasker.tbl" from repeatmaker_table
		file "ltr.gff3" from ltr_gff 
	output:
		file "Repeat.gff3" into RepeatGFFResult
		file "Scaffolds.Masked.fasta" into RepeatSeqResult
		file "*.xls" into ResultTable
	
	
	"""
		 cat Denovo.gff3  Protein.gff3    ltr.gff3> combined.gff
		 bedtools  sort  -i combined.gff  >All.gff
		 bedtools  merge -c 2,3,6,7,8,9  -o distinct -i All.gff  | cut -f 1,4,5,2,3,6,7,8,9 >merged.gff
		 
			 
		cp RepeatMasker.tbl RepeatMasker.xls
		RepeatRename.sh merged.gff > Repeat.gff3 
		RepeatMaskerSequenceFromGFF.py -i Scaffolds.fa  -g Repeat.gff3  -o Scaffolds.Masked.fasta
	"""
}
 


 
 
 
 
 
 
 
 
 
 
 
 
 
