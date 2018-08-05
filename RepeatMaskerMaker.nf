#!/usr/bin/env nextflow
VERSION = "1.0.0"
/*
#params.trnaprot = 'http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz'
#params.trnanuc = '/home/nfs/Database/tRNA/eukaryotic-tRNAs.fa'
params.path = '/home/nfs/SOFTWARE/Other/nf-repeatmasking/bin'
*/
params.reference = 'data/example/genome.fasta'
params.trnaprot = '/home/nfs/Database/tRNA/Tpases020812.gz'
params.trnanuc = '/home/nfs/Database/tRNA/eukaryotic-tRNAs.fa.gz'
params.path = '/home/nfs/SOFTWARE/Other/nf-repeatmasking/bin'
params.outdir = 'output'
params.minRepLength = 250

trnanuc = file(params.trnanuc)
trnaprot = file(params.trnaprot)
reference = file(params.reference)

log.info ""
log.info "R E P E A T M A S K I N G  ~  version " + VERSION
log.info "reference input    : ${params.reference}"
log.info "output directory   : ${params.outdir}"
log.info ""
process MITE_Hunter{
 executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"
  
  input:
	file 'genome.fasta' from reference
  output:
	file "*.lib" into MITELib
	"""
	MITE_Hunter_manager.pl -i genome.fasta -g MITE -n 32 -S 12345678
cat *Step8_*.fa > MITE.lib
	"""
}

process Index{
 executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"
  
  input:
	file 'genome.fasta' from reference
  output:
	file("CTG.*") into Index
	"""
	gt suffixerator -db genome.fasta -indexname CTG -tis -suf -lcp -des -ssp -dna
	"""
}



process recentLTRs {
  executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"
  cache 'deep'

  input:
  file 'genome.fasta' from reference
  file 'eukaryotic-tRNAs.fasta.gz' from trnanuc
  file ctg_index from Index 

  output:
  set age, 'seqfile.result' into ltrHarvestResultsNew
  set age, 'seqfile.result' into ltrHarvestResultsForExemplarNew
  set age, 'seqfile.outinner' into ltrHarvestInnerNew
  set age, 'seqfile.outinner' into outinnerForBlastXNew
  set age, 'CRL_Step2_Passed_Elements.fasta', 'Repeat_down*.fasta', 'Repeat_up*.fasta' into recentLTRs

  script:
  age = 'new'
  """


	gt ltrharvest \
	 -index CTG \
	 -out seqfile.out \
	 -outinner seqfile.outinner \
	 -gff3 seqfile.gff \
	 -minlenltr 100 \
	 -maxlenltr 6000 \
	 -mindistltr 1500 \
	 -maxdistltr 25000 \
	 -mintsd 5 \
	 -maxtsd 5 \
	 -motif tgca \
	 -similar 99 \
	 -vic 10 \
	> seqfile.result

	gt gff3 \
	 -sort seqfile.gff \
	> seqfile.gff.sort

	zcat eukaryotic-tRNAs.fasta.gz > eukaryotic-tRNAs.fasta

	gt ltrdigest \
	 -trnas eukaryotic-tRNAs.fasta \
	 seqfile.gff.sort \
	 CTG \
	> seqfile.gff.dgt

	CRL_Step1.pl \
	 --gff seqfile.gff.dgt

	CRL_Step2.pl \
	 --step1 CRL_Step1_Passed_Elements.txt \
	 --repeatfile seqfile.out \
	 --resultfile seqfile.result \
	 --sequencefile genome.fasta \
	 --removed_repeats CRL_Step2_Passed_Elements.fasta
  """
}

process olderLTRs {
  executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  cache 'deep'

  input:
  file ctg_index from Index 
  file 'genome.fasta' from reference
  file 'eukaryotic-tRNAs.fasta.gz' from trnanuc

  output:
  set age, 'seqfile.result' into ltrHarvestResultsOld
  set age, 'seqfile.result' into ltrHarvestResultsForExemplarOld
  set age, 'seqfile.outinner' into ltrHarvestInnerOld
  set age, 'seqfile.outinner' into outinnerForBlastXOld
  set age, 'CRL_Step2_Passed_Elements.fasta', 'Repeat_down*.fasta', 'Repeat_up*.fasta' into olderLTRs

  script:
  age = 'old'
  """


gt ltrharvest \
 -index CTG \
 -out seqfile.out \
 -outinner seqfile.outinner \
 -gff3 seqfile.gff \
 -minlenltr 100 \
 -maxlenltr 6000 \
 -mindistltr 1500 \
 -maxdistltr 25000 \
 -mintsd 5 \
 -maxtsd 5 \
 -vic 10 \
> seqfile.result

gt gff3 \
 -sort seqfile.gff \
> seqfile.gff.sort

zcat eukaryotic-tRNAs.fasta.gz > eukaryotic-tRNAs.fasta

gt ltrdigest \
 -trnas eukaryotic-tRNAs.fasta \
 seqfile.gff.sort \
 CTG \
> seqfile.gff.dgt

CRL_Step1.pl \
 --gff seqfile.gff.dgt

CRL_Step2.pl \
 --step1 CRL_Step1_Passed_Elements.txt \
 --repeatfile seqfile.out \
 --resultfile seqfile.result \
 --sequencefile genome.fasta \
 --removed_repeats CRL_Step2_Passed_Elements.fasta
  """
}

ltrs = recentLTRs.mix(olderLTRs)
ltrHarvestResults = ltrHarvestResultsOld.mix(ltrHarvestResultsNew)
ltrHarvestInner = ltrHarvestInnerOld.mix(ltrHarvestInnerNew)
outinnerForBlastX = outinnerForBlastXOld.mix(outinnerForBlastXNew)
ltrHarvestResultsForExemplar = ltrHarvestResultsForExemplarOld.mix(ltrHarvestResultsForExemplarNew)

process CRL_Step3 {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  tag { age }

  input:
  set age, 'CRL_Step2_Passed_Elements.fasta', 'Repeat_down*.fasta', 'Repeat_up*.fasta' from ltrs

  output:
  set age, 'CRL_Step3_Passed_Elements.fasta' into step3Passed
  set age, 'CRL_Step3_Passed_Elements.fasta' into step3PassedForExemplars

  """
CRL_Step3.pl \
 --directory . \
 --step2 CRL_Step2_Passed_Elements.fasta \
 --pidentity 60 \
 --seq_c 25
  """
}

ltrHarvestResults
.combine(step3Passed, by: 0)
.set { nestedInput }

process identifyNestedInsetions {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  tag { age }
  input:
  file 'genome.fasta' from reference
  set age, 'seqfile.result', 'CRL_Step3_Passed_Elements.fasta' from nestedInput
  file "MITE.lib" from MITELib
  output:
  set age, 'repeats_to_mask_LTR.fasta' into repeatsToMaskLTR

  """
ltr_library.pl \
 --resultfile seqfile.result \
 --step3 CRL_Step3_Passed_Elements.fasta \
 --sequencefile genome.fasta
cat lLTR_Only.lib \
| awk 'BEGIN {RS = ">" ; FS = "\\n" ; ORS = ""} \$2 {print ">"\$0}' \
> repeats_to_mask_LTR.fasta

cat MITE.lib>> repeats_to_mask_LTR.fasta
  """
}

process RepeatMasker1 {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  tag { age }

  input:
  set age, 'repeats_to_mask_LTR.fasta', 'seqfile.outinner' from repeatsToMaskLTR.combine(ltrHarvestInner, by: 0)

  output:
  set age, 'seqfile.outinner.out', 'seqfile.outinner.masked' into repeatMasker1Unclean

  """
RepeatMasker \
 -lib repeats_to_mask_LTR.fasta \
 -nolow \
 -no_is \
 -dir . \
 seqfile.outinner

if [ ! -f seqfile.outinner.masked ]; then
  cp seqfile.outinner seqfile.outinner.masked
fi
  """
}

process cleanRM {
  tag { age }
  executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  input:
  set age, 'seqfile.outinner.out', 'seqfile.outinner.masked' from repeatMasker1Unclean

  output:
  set age, 'seqfile.outinner.clean' into repeatMasker1Clean

  """
cleanRM.pl seqfile.outinner.out seqfile.outinner.masked > seqfile.outinner.unmasked
rmshortinner.pl seqfile.outinner.unmasked 50 > seqfile.outinner.clean
  """
}

process blastX {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  tag { age }
  cpus 32

  input:
  file 'Tpases020812DNA.fasta.gz' from trnaprot
  set age, 'seqfile.outinner.clean', 'seqfile.outinner' from repeatMasker1Clean.combine(outinnerForBlastX, by: 0)

  output:
  set age, 'passed_outinner_sequence.fasta' into blastxPassed

  """
zcat Tpases020812DNA.fasta.gz > Tpases020812DNA.fasta
makeblastdb -in Tpases020812DNA.fasta -dbtype prot
blastx \
 -query seqfile.outinner.clean \
 -db Tpases020812DNA.fasta \
 -evalue 1e-20 \
 -num_descriptions 10 \
 -num_threads ${task.cpus} \
 -out seqfile.outinner.clean_blastx.out.txt

outinner_blastx_parse.pl \
 --blastx seqfile.outinner.clean_blastx.out.txt \
 --outinner seqfile.outinner

if [ ! -s passed_outinner_sequence.fasta ]; then
  echo -e '>dummy empty sequence\nACTACTAC' > passed_outinner_sequence.fasta
fi
  """
}

blastxPassed
.combine(step3PassedForExemplars, by: 0)
.combine(ltrHarvestResultsForExemplar, by: 0)
.set { forExemplarBuilding }

process buildExemplars {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  tag { age }
  cpus 16

  input:
  file 'genome.fasta' from reference
  set age, 'passed_outinner_sequence.fasta', 'CRL_Step3_Passed_Elements.fasta', 'seqfile.result' from forExemplarBuilding

  output:
  set age, 'LTR.lib' into exemplars

  """
CRL_Step4.pl \
 --step3 CRL_Step3_Passed_Elements.fasta \
 --resultfile seqfile.result \
 --innerfile passed_outinner_sequence.fasta \
 --sequencefile genome.fasta

for lib in lLTRs_Seq_For_BLAST.fasta Inner_Seq_For_BLAST.fasta; do
  makeblastdb -in \$lib -dbtype nucl
  blastn \
   -query \${lib} \
   -db \${lib} \
   -evalue 1e-10 \
   -num_threads ${task.cpus} \
   -num_descriptions 1000 \
   -out \${lib}.out
done

CRL_Step5.pl \
 --LTR_blast lLTRs_Seq_For_BLAST.fasta.out \
 --inner_blast Inner_Seq_For_BLAST.fasta.out \
 --step3 CRL_Step3_Passed_Elements.fasta \
 --final LTR.lib \
 --pcoverage 90 \
 --pidentity 80
  """
}

newLTRs = Channel.create()
oldLTRs = Channel.create()

exemplars
.choice( newLTRs, oldLTRs ) { it[0] == "new" ? 0 : 1 }

process removeDuplicates {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"


  input:
  set _, 'ltrs.new.fasta' from newLTRs
  set _, 'ltrs.old.fasta' from oldLTRs

  output:
  set 'ltrs.old.fasta.masked', 'ltrs.new.fasta' into bothLTRforMasking
  """RepeatMasker -lib ltrs.new.fasta -dir . ltrs.old.fasta
  if [ ! -f ltrs.old.fasta.masked ]; then
  cp ltrs.old.fasta ltrs.old.fasta.masked
fi"""
  
}

process filterOldLTRs {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"


  input:
  set 'ltrs.old.fasta.masked', 'ltrs.new.fasta' from bothLTRforMasking

  output:
  file 'allLTRs.fasta' into allLTR

  """
remove_masked_sequence.pl \
--masked_elements ltrs.old.fasta.masked \
--outfile ltrs.old.final.fasta
cat ltrs.new.fasta ltrs.old.final.fasta > allLTRs.fasta
  """
}

allLTR
.splitFasta(record: [id: true, sequence: true ])
.collectFile( name: 'allLTRs.fasta' ) { ">" + it.id + "#LTR\n" + it.sequence }
.tap { allLTR2 }
.set { inputForRM2 }

process RepeatMasker2 {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  cpus 10

  input:
  file 'genome.fasta' from reference
  file 'allLTR.lib' from inputForRM2

  output:
  file 'genome.fasta.masked' into genomeLtrMasked

  """
RepeatMasker \
 -no_is \
 -nolow \
 -pa ${task.cpus} \
 -lib allLTR.lib \
 -dir . \
 genome.fasta
  """
}

process RepeatModeler {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  cpus 10

  input:
  file 'genome.masked' from genomeLtrMasked

  output:
  file 'consensi.fa.classified' into rmOutput

  """
rmaskedpart.pl genome.masked 50 | sed '/^\$/d' > umseqfile
BuildDatabase -name umseqfiledb -engine ncbi umseqfile
RepeatModeler -pa ${task.cpus} -database umseqfiledb >& umseqfile.out
mv RM*/consensi.fa.classified consensi.fa.classified
  """
}

identityUnknown = Channel.create()
identityKnown = Channel.create()

rmOutput
.splitFasta(record: [id: true, text: true])
.choice(identityUnknown, identityKnown) { record -> record.id =~ /#Unknown/ ? 0 : 1 }

identityUnknown
.collectFile() { record -> ['unknown.fasta', record.text] }
.set { repeatmaskerUnknowns }

identityKnown
.collectFile() { record -> ['known.fasta', record.text] }
.tap { repeatmaskerKnowns1 }
.set{ repeatmaskerKnowns }

process searchForUnidentifiedElements {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"


  input:
  file 'genome.fasta' from reference
  file 'unknown_elements.fasta' from repeatmaskerKnowns1

  output:
  set 'genome.fasta.align', 'unknown_elements.fasta' into unknownAlignments

  """
RepeatMasker \
 -lib unknown_elements.fasta \
 -alignments \
 -nolow \
 -no_is \
 -dir . \
 -inv \
 genome.fasta
  """
}

process derip {
    executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=1 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"


  input:
  set 'genome.fasta.align', 'unknown_elements.fasta' from unknownAlignments

  output:
  file 'deripped.unknowns.fasta' into derippedUnknowns

  """
rmalignment_to_fasta.rb genome.fasta.align unknown_elements.fasta
for file in alignment*; do
  derip.rb \$file
done > deripped.unknowns.fasta
  """
}

process classifyDeripped {
  executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"


  input:
  file 'transposases.fasta.gz' from trnaprot
  file 'repeatmodeler_unknowns.deripped.fasta' from derippedUnknowns

  output:
  file 'identified_elements.txt' into identifiedDerippedTransposons
  file 'unknown_elements.txt' into unidentifiedDerippedTransposons

  """
zcat transposases.fasta.gz > transposases.fasta
makeblastdb \
 -in transposases.fasta \
 -dbtype prot
blastx \
 -query repeatmodeler_unknowns.deripped.fasta \
 -db transposases.fasta \
 -evalue 1e-10 \
 -num_descriptions 10 \
 -num_threads ${task.cpus} \
 -out modelerunknown_blast_results.txt
transposon_blast_parse.pl \
 --blastx modelerunknown_blast_results.txt \
 --modelerunknown repeatmodeler_unknowns.deripped.fasta
  """
}

identifiedDerippedTransposons.subscribe{ println("Identified, deripped: ${it}") }

process transposonBlast {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  cpus 4

  input:
  file 'transposases.fasta.gz' from trnaprot
  file 'repeatmodeler_unknowns.fasta' from repeatmaskerUnknowns

  output:
  file 'identified_elements.txt' into identifiedTransposons
  file 'unknown_elements.txt' into unidentifiedTransposons

  """
zcat transposases.fasta.gz > transposases.fasta
makeblastdb \
 -in transposases.fasta \
 -dbtype prot
blastx \
 -query repeatmodeler_unknowns.fasta \
 -db transposases.fasta \
 -evalue 1e-10 \
 -num_descriptions 10 \
 -num_threads ${task.cpus} \
 -out modelerunknown_blast_results.txt
transposon_blast_parse.pl \
 --blastx modelerunknown_blast_results.txt \
 --modelerunknown repeatmodeler_unknowns.fasta
  """
}

process predict_ncRNA {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=32 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"


  input:
  file 'genome.fasta' from reference

  output:
  file 'cm_out' into cmtblouts

  """
cp /home/nfs/Database/rnas.cm models.cm
cmpress -F models.cm
cmsearch --cpu 32 --tblout cm_out --cut_ga models.cm genome.fasta
  """
}

process merge_ncrnas {
  publishDir "${params.outdir}/noncoding-rna", mode: 'copy'
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  cache 'deep'

  input:
  file 'cmtblout' from cmtblouts

  output:
  file 'ncrna.gff3' into ncrnafile

  """
  infernal_to_gff3.lua < ${cmtblout} > 1
  gt gff3 -sort -tidy -retainids 1 > ncrna.gff3
  """
}

ncrnafile.println()

repeatmaskerKnowns
.mix(identifiedTransposons)
.collectFile() { it.text }
.combine(allLTR2)
.set { knownRepeats }

process repeatMaskerKnowns {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  publishDir "${params.outdir}/repeatMaskerKnowns", mode: 'copy'

  input:
  file 'reference.fasta' from reference
  set 'knownTransposons.lib', 'allLTRs.lib' from knownRepeats

  output:
  set 'reference.fasta.out', 'reference.fasta.masked' into repeatMaskerKnownsMasked
  set 'reference.fasta.out.gff', 'knownRepeats.library.fasta'

  """
cat *.lib > knownRepeats.library.fasta
RepeatMasker \
 -lib knownRepeats.library.fasta \
 -nolow \
 -no_is \
 -dir . \
 -gff \
 -s \
 reference.fasta
  """
}

process removeShortMatches {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"

  publishDir "${params.outdir}/cleanMasked", mode: 'copy'

  input:
  file 'reference.fa' from reference
  set 'rm.out', 'rm.masked' from repeatMaskerKnownsMasked

  output:
  set 'rm.trimmed.out', 'rm.trimmed.masked' into repeatMaskerKnownsMaskedTrimmed

  """
head -n 3 rm.out > rm.trimmed.out
tail -n +4 rm.out | awk '\$7 - \$6 > ${params.minRepLength}' >> rm.trimmed.out
tail -n +4 rm.out | awk 'BEGIN{OFS="\\t"} \$7 - \$6 > ${params.minRepLength} {print \$5, \$6, \$7}' >> rm.trimmed.bed
maskFastaFromBed -fi reference.fa -bed rm.trimmed.bed -fo rm.trimmed.masked -soft
  """
}

process octfta {
   executor 'pbs'

  clusterOptions  " -d $PWD  -l nodes=1:ppn=16 -v PATH=$PATH:${params.path};PERLLIB=$PERLLIB"


  input:
  file 'reference.fa' from reference
  set 'rm.out', 'rm.masked' from repeatMaskerKnownsMaskedTrimmed

  output:
  file 'summary.tsv' into repeatmaskerSummaryTable

  """
build_dictionary.pl --rm rm.out > ltr.dict
one_code_to_find_them_all.pl --rm rm.out --ltr ltr.dict --fasta reference.fa
echo -e 'Family\\tElement Length\\tFragments\\tCopies\\tSolo_LTR\\tTotal_Bp\\tCover\\tchrname' > summary.tsv
for file in *.copynumber.csv; do
  chrname=`echo \$file | sed -e 's/^rm\\.out_//' -e 's/.copynumber.csv\$//'`
  awk -v chrname=\$chrname 'BEGIN{OFS="\\t"} NR>1 && /^[^#]/ {print(\$0, chrname)}' \$file
done >> summary.tsv
  """
}

process summarise {
  publishDir "${params.outdir}/summarise", mode: 'copy'

  input:
  file 'summary.tsv' from repeatmaskerSummaryTable

  output:
  set 'summary.bycontig.tidy.tsv', 'summary.tidy.tsv'

  """
#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)

data <- read.table('summary.tsv', header=TRUE) %>%
        separate(Family, into=c("Family", "Subfamily"), sep="/") %>%
        group_by(chrname, Family, Subfamily) %>%
        summarise(fragment.count = sum(Fragments), length = sum(Total_Bp)) %>%
        unite("Family", Family, Subfamily, sep="/")

write.table(data, file='summary.bycontig.tidy.tsv')

data <- read.table('summary.tsv', header=TRUE) %>%
        separate(Family, into=c("Family", "Subfamily"), sep="/") %>%
        group_by(Family, Subfamily) %>%
        summarise(fragment.count = sum(Fragments), length = sum(Total_Bp)) %>%
        unite("Family", Family, Subfamily, sep="/")

write.table(data, file='summary.tidy.tsv', row.names = FALSE)
  """
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
