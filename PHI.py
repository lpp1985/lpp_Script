#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/11/24
"""

from lpp import * 
from optparse import OptionParser
from termcolor import colored
import os


if __name__ == '__main__':

	usage = '''usage: python %prog [options] 

	It can automaticly run BWA!!'''
	parser = OptionParser(usage =usage ) 

	parser.add_option("-i", "--Input", action="store", 
	                  dest= "Input",
	                  type='string',
	                  help="the input seq")

	parser.add_option("-o", "--Output", action="store", 
	                  dest="output",
	                  type='string',
	                  help="output")
	parser.add_option("-x", "--XLS", action="store", 
	                  dest="xls",
	                  type='string',
	                  help="vfdb xls")	

	(options, args) = parser.parse_args()
	RAW = open( options.xls, 'rU')
	all_info = {}
	for line in RAW:
		line_l = line.split("\n")[0].split("\t")
		all_info[ line_l[3]] = line
	Input = options.Input
	OUTPUT = open(options.output, 'w')
	OUTPUT.write( "Name\tPHI_Annotation\tPHI_Identity\tPHI_Evalue\tCuration_comments_temp\t_to_do_\tRecord ID	PHI_MolConn_ID	IdentifierTypeOfProteinID	ProteinID	IdentifierTypeOfGeneLocusID	GeneLocusID	AA sequence #no EMBL#	NT sequence #no EMBL#	Genomic sequence providing strain	Gene_name	Genome location	Specific modification/s to the targeted protein or promoter	Accession ID for the modified genetic element	Known interacting protein(s) in the pathogen	Interacting protein - locus ID	Multiple_mutation	Pathogen_NCBI_species_Taxonomy ID	Pathogen_species	Pathogen_NCBI_strain_Taxonomy_ID	Experimental_strain	Disease_name	Host_descripton	Host_NCBI_Taxonomy_ID	Experimental_host_species	Host_strain_genotype_or cultivar_taxonomy_ID_NCBI	HostGenotype_definedGene of interest	HostGenotype_definedGene of interest _AccessionID	tissue_type	Function	GO_annotation	Database	Pathway_secretion_systems	Phenotype_of_mutant	Mating_defect_prior_to_penetration	Pre_penetration_defect	Penetration_defect	Post_penetration_defect	Disease_development_macroscopically_visible	Vegetative_spores	Sexual_spores	In_vitro_growth	Spore_germination	Essential_gene_Lethal_knockout	Inducer	ChemicalAccession (Chebi/CAS)	Tested Host_target	TestedHostTarget_AccessionID	Interaction phenotype	Host_response	Experimental_evidence	Transient Assay Experimental Evidence	Species_Expert	Entered_by	Literature_ID	Literature_source	DOI	Full_citation	Author_email	Comments	Reference 	Year_published	Curation details	File_name_pdf_files_provided	_batch_no_	Curation date	_curator_organisation_	__lab__	__FG_mycotoxin__	Additional IdentifierTypeOfGeneLocusID	Additional GeneLocusID	Anti-infective (Chemical)	Compound	Target site	Group name	Chemical group	Mode in planta	Mode of action 	FRAC CODE	Additional comments  on anti-infectives	Cellular/sub-cellular phenotypes reported in paper 'yes' or 'no'	Horizontal gene transfer 'yes' or 'no'	phi3_pars	AC For First Host Targets. Pathogen gene maps in EnsemblFungi/protist?	AC For First Host Targets. Host gene maps in EnsemblPlants?\n")
	all_result = os.popen( """ blastp  -query %s -db /home/nfs/Database/PHI/PHI  -evalue 1e-35 -outfmt  "6 qseqid  stitle pident  evalue" -num_threads 16  -num_alignments  1  """ % (Input))
	for line in all_result:
		data = re.search("(PHI\:\d+)", line).group(1)
		OUTPUT.write(line.strip() + '\t' + all_info[data] + '\n')
