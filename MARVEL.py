#!/usr/bin/env python

import multiprocessing
import marvel
import marvel.config
import marvel.queue
import sys
import os,re
### settings

DB         = sys.argv[1]
COVERAGE   = sys.argv[2]
PATH="{path}/"
PATH_SCRIPTS="{path_scripts}"
DB_FIX     = DB + "_FIX"
mem = os.popen( "grep \"MemTotal\" /proc/meminfo " )
mem = int(re.search("(\d+)",mem.read()).group(1))
mem = mem/10/1024/1024

PARALLEL   = mem

### patch raw reads

q = marvel.queue.queue(DB, COVERAGE, PARALLEL)

### run daligner to create initial overlaps
q.plan("{db}.dalign.plan")

### run LAmerge to merge overlap blocks
q.plan("{db}.merge.plan")

# create quality and trim annotation (tracks) for each overlap block
q.block("{path}//LAq -b {block} {db} {db}.{block}.las")
# merge quality and trim tracks
q.single("{path}//TKmerge -d {db} q")
q.single("{path}//TKmerge -d {db} trim")

# run LAfix to patch reads based on overlaps
q.block("{path}//LAfix -g -1 {db} {db}.{block}.las {db}.{block}.fixed.fasta")
# join all fixed fasta files
q.single("!cat {db}.*.fixed.fasta > {db}.fixed.fasta")

# create a new Database of fixed reads (-j numOfThreads, -g genome size)
q.single("{path_scripts}/DBprepare.py -s 50 -r 2 -j 4 -g 4600000 {db_fixed} {db}.fixed.fasta", db_fixed = DB_FIX)

q.process()

##### assemble patched reads

q = marvel.queue.queue(DB_FIX, COVERAGE, PARALLEL)

### run daligner to create overlaps
q.plan("{db}.dalign.plan")
### run LAmerge to merge overlap blocks
q.plan("{db}.merge.plan")

### start scrubbing pipeline

##### for larger genomes (> 100MB) LAstitch can be run with the -L option (preload reads)
##### with the -L option two passes over the overlap files are performed:
##### first to buffer all reads and a second time to stitch them
##### otherwise the random file access can make LAstitch pretty slow.
##### Another option would be, to keep the whole db in cache (/dev/shm/)
q.block("{path}//LAstitch -L -f 50 {db} {db}.{block}.las {db}.{block}.stitch.las")

##### create quality and trim annotation (tracks) for each overlap block and merge them
q.block("{path}//LAq -s 5 -T trim0 -b {block} {db} {db}.{block}.stitch.las")
q.single("{path}//TKmerge -d {db} q")
q.single("{path}//TKmerge -d {db} trim0")

##### create a repeat annotation (tracks) for each overlap block and merge them

q.block("{path}/LAfilter -o 2000 -L -p -t trim0   {db} {db}.{block}.stitch.las {db}.{block}.ol.las")
q.block("{path}/LArepeat -l 1.7 -h  2.0  -t repeats -c {coverage} -b {block} {db} {db}.{block}.ol.las")
q.single("{path}/TKmerge -d {db} repeats")
q.block("  {path}/TKhomogenize -i repeats -I hrepeats -b {block} {db} {db}.{block}.ol.las ")
q.single(" {path}/TKcombine -d {db} hrepeats \#.hrepeats")
q.single(" {path}/TKcombine {db} frepeats repeats hrepeats")
##### detects "borders" in overlaps due to bad regions within reads that were not detected
##### in LAfix. Those regions can be quite big (several Kb). If gaps within a read are
##### detected LAgap chooses the longer part oft the read as valid range. The other part(s) are
##### discarded
##### option -L (see LAstitch) is also available
q.block("{path}//LAgap -L -t trim0 {db} {db}.{block}.ol.las {db}.{block}.gap.las")

##### create a new trim1 track, (trim0 is kept)
q.block("{path}//LAq -s 5 -u -t trim0 -T trim1 -b {block} {db} {db}.{block}.gap.las")
q.single("{path}//TKmerge -d {db} trim1")

##### based on different filter critera filter out: local-alignments, repeat induced-alifnments
##### previously discarded alignments, ....
##### -r repeats, -t trim1 ... use repeats and trim1 track
##### -n 500  ... overlaps must be anchored by at least 500 bases (non-repeats)
##### -u 0    ... overlaps with unaligned bases according to the trim1 interval are discarded
##### -o 2000 ... overlaps shorter than 2k bases are discarded
##### -p      ... purge overlaps, overlaps are not written into the output file
##### option -L (see LAstitch) is also available
q.block("{path}//LAfilter -L -n 100 -r repeats -t trim1 -T -o 2000 -u 0 {db} {db}.{block}.gap.las {db}.{block}.filtered.las")

##### merge all filtered overlap files into one overlap file
q.single("{path}//LAmerge -S filtered {db} {db}.filtered.las")

##### create overlap graph
q.single("{path}//OGbuild -t trim1 {db} {db}.filtered.las {db}.graphml")

##### tour the overlap graph and create contigs paths
q.single("{path_scripts}/OGtour.py -c {db} {db}.graphml")

q.single("{path}//LAcorrect -j 4 -r {db}.tour.rids {db} {db}.filtered.las {db}.corrected")
q.single( "cat {db}.corrected.*.fasta >Corrected.fa" )
q.single("{path}//FA2db -c {db}_CORRECTED Corrected.fa")

##### create contig fasta files
q.single("{path_scripts}/tour2fasta.py -c {db}_CORRECTED -t trim1 {db} {db}.tour.graphml {db}.tour.paths")

### optional: create a layout of the overlap graph which can viewed in a Browser (svg) or Gephi (dot)
q.single("{path}//OGlayout -R {db}.tour.graphml {db}.tour.layout.svg")
q.single("{path}//OGlayout -R {db}.tour.graphml {db}.tour.layout.dot")

q.process()
