#!/home/nfs/SOFTWARE/bin/nextflow
params.Nucl = "nucl.fa"
params.Protein = "Protein.fa"
params.out = "./Annotation"
params.eval = "1e-5"
name_prefix = file(params.Protein).getBaseName()
output_prefix= file(params.out+'/Detail/'+name_prefix).toAbsolutePath()

Channel.fromPath(params.Protein).into { pro_kegg; pro_go;pro_cog; pro_swiss; pro_nr }
evalue = params.eval 

Channel.fromPath(params.Nucl).into { nucl_kegg; nucl_nt }

process Nr_Annotation {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
		file "pro.pep" from pro_nr

 
    """
	 Nr_Annotation.py -i pro.pep -o ${output_prefix}/Nr/${name_prefix} -e ${evalue}
    
    """
}

process Nt_Annotation {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
		file "query.fa" from nucl_nt

 
    """
	 Nt_Annotation.py -i query.fa   -o ${output_prefix}/Nt/${name_prefix}.xls -e ${evalue}

    
    """
} 
 
process Swiss_Annotation {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
		file "pro.pep" from pro_swiss

 
    """
	 Swiss_Annotation.py -i pro.pep -o ${output_prefix}/Swiss/${name_prefix} -e ${evalue}
    
    """
}


process GO_Annotation {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
		file "pro.pep" from pro_go

 
    """
	 GO_Annotation.py -i pro.pep -o ${output_prefix}/GO/${name_prefix} -e ${evalue}
    
    """
}

 
 
process COG_Annotation {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
		file "pro.pep" from pro_cog

 
    """
	 COG_Annotation.py -i pro.pep -o ${output_prefix}/eggNOG/${name_prefix} -e ${evalue}
   
    """
}
 
 
process KEGG_Annotation {
    executor 'pbs'

    clusterOptions  " -d $PWD  -l nodes=1:ppn=8 -v PATH=$PATH"
    input:
		file "pro.pep" from pro_kegg
		file "nucl.fa" from nucl_kegg

 
    """
	KEGG_Annotation.py -p pro.pep  -n nucl.fa  -o  ${output_prefix}/KEGG/${name_prefix}  -e ${evalue}

   
    """
}
 
 
 
 
 
