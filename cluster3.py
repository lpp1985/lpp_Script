#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/10/11
"""
import networkx as nx
from lpp import *
RAW = open(sys.argv[1],'rU')
graph = nx.Graph()
FASTA = fasta_check( open(sys.argv[2],'rU') )
length_node = {}
seq_hash = {}
for t,s in FASTA:
  name = t[1:].split()[0]
  graph.add_node(name)
  length_node[name] = len(s)
  seq_hash[name] = s
  
for line in RAW:
  line_l = line.strip().split("\t")
  query,subject = line_l[11],line_l[12]
  if "[C" in line or "[I" in line  and query !=subject:
    graph.add_edge(query,subject)
  query_perc = float(line_l[6])*max([ float(line_l[9]),  float(line_l[10])  ])
  if query_perc>9000:
    graph.add_edge(query,subject)

END = open("cluster_out.clust",'w')
FASTA = open("cluster_repre.fasta",'w')
LIST = open("cluster_seq.list",'w')
i=0
def max_length(data_list):
  candidate_list = sorted(data_list,key=lambda x:length_node[x] )[::-1]
  candidate = candidate_list[0]
  if filter(lambda x: 'round' in x,data_list):
    for i in candidate_list:
      if 'round' in i:
        candidate=i
        break
      
  
  sequ = seq_hash[candidate]
  return candidate,sequ
for key in nx.connected_components(graph):
  i+=1
  END.write("Cluster%s\t%s\n"%(i,'; '.join(key)))
  candidate,sequ = max_length(key)
  FASTA.write("Cluster%s\n%s"%(i,sequ))
  LIST.write("Cluster%s\t%s\n"%(i,candidate))

  