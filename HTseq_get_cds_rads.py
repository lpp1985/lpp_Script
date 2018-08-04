#!/usr/bin/python
import HTSeq  ,sys

BAM = HTSeq.BAM_Reader( sys.argv[1]   )
GTF = HTSeq.GFF_Reader( sys.argv[2],end_included=True )
#BAM_WRITER = HTSeq.BAM_Writer.from_BAM_Reader( sys.argv[3], BAM )
cds_all = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
counts = {}
all_cds = {}
for feature in GTF:

    if feature.type=='exon':
	#print( feature.attr )
	name = feature.attr[ 'transcript_id' ]
	cds_all[ feature.iv   ] += name
	counts[ name ] = 0
    if feature.type =='CDS':
	all_cds[  name] = ''

#print( counts )
#for a,b in cds_all.steps():
#    print(a,b)
#print list( cds_all[  HTSeq.GenomicInterval('MG1655', 14436,14531 ) ]  .steps() )
for alnmt in BAM:
    if alnmt.aligned and alnmt.proper_pair:
	#intersection_set = None
	#for iv2, step_set in cds_all[ alnmt.iv ].steps():
	    #if intersection_set is None:
		#intersection_set = step_set.copy()		
	    #else:
		
		#intersection_set.intersection_update( step_set )
	#print step_set
	#print( cds_all[ alnmt.iv ].steps().next() )
	gene_set = cds_all[ alnmt.iv ].steps().next()[-1]
	if len( gene_set ) == 1:
	#    print( gene_set )
	    if  list(gene_set)[0] in all_cds: 

		counts[ list(gene_set)[0] ] += 1
#		BAM_WRITER.write( alnmt)
#BAM_WRITER.close()
for gene,value in counts.items():
    print( '%s\t%s'%( gene,value  ) )
