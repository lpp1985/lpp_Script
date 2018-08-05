#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/29
"""
from lpp import *
import networkx as nx
from lpp import *
import string
libary=string.maketrans('+-','-+')
Graph = Ddict()
def reverse(data):
    #转换成为反向节点
    return [x.translate(libary) for x in data[::-1] ]

def add_bi_edge(start,end,weight):
    global Graph
    
    Graph[start][end]["weight"] = weight
    start,end = reverse( [start,end] )
    Graph[start][end]["weight"] = weight
    
def remove_bi_edge(start,end):
    if start in Graph and end in Graph[start]:
        del Graph[start][end]
    start,end = reverse( [start,end] )
    if start in Graph and end in Graph[start]:
        del Graph[start][end]    




END = open("true5.ovl",'w')
END2 = open("filter5.ovl",'w')


contain={}
#for l in open(sys.argv[1]):
    #line = l.strip()
    #line_l = line.split()
    #if line.endswith("contained"):
    #contain[line_l[0]] = ""
    #END2.write(l)
    #elif line.endswith("contains"):
    #contain[line_l[1]] = ""
    #END2.write(l)
add_edge = {}
for line in open(sys.argv[1]):

    line = line.strip()
    if line.endswith("overlap"):
        line_l = line.split()
        if line_l[0] in contain or line_l[1] in contain:
            continue
        score = int(line_l[2])
        query_start = int(line_l[5])
        query_end = int(line_l[6])
        query_length = int(line_l[7])
        frame =  int(line_l[8])
        sub_start = int(line_l[9])
        sub_end = int(line_l[10])
        sub_length = int(line_l[11])
        query_name = line_l[0]+"+"
        subject_name= line_l[1]
        if line_l[0] == line_l[1]:
            continue
        if frame==1:
            dir="-"
        else:
            dir="+"
        subject_name +=dir        
        if int(line_l[5])<25:
            if frame==0:
                if sub_length - sub_end<50:
                    END.write(line+'\n')
                    add_edge[ (subject_name,query_name, score)   ] = ""
                    #Graph.add_bi_edge( subject_name,query_name, score)
            elif frame==1:
                if sub_end<50:
                    END.write(line+'\n')
                    add_edge[ (subject_name,query_name, score)   ] = ""
                    #Graph.add_bi_edge( subject_name,query_name, score)
        elif query_length-query_end<50:
            if frame==0:
                if sub_end<50:
                    END.write(line+'\n')	
                    add_edge[ (query_name,subject_name, score)   ] = ""
                    #Graph.add_bi_edge( query_name,subject_name, score)
            elif frame==1:
                if sub_length - sub_end<50:
                    END.write(line+'\n')
                    add_edge[ (query_name,subject_name, score)   ] = ""
                    #Graph.add_bi_edge(query_name,subject_name,  score)
i=0
for q,s,sc in add_edge:

    add_bi_edge(q,s,  sc)
END.close()
remove_list = []

#Delete Orphan Node if other candidate exists

for node in Graph:
    cache_list = []
    tag = 0
    for each_pred in Graph[node]:
        if each_pred not in Graph:
            cache_list.append( (node,each_pred) )
            continue
        for pre_pred in Graph[each_pred]:
            if len(pre_pred)>0:
                tag = 1
    if tag==1:
        remove_list+=cache_list
for s,t in set(remove_list):
    remove_bi_edge(s,t)
#每个节点保留两个最好的overlap
New_Graph = Ddict()
for node in Graph:
    best_two = sorted( Graph[query_name],key = lambda x: Graph[query_name][x]["weight"] )[:2]
    if len(best_two)>1:
        
        if float(Graph[query_name][best_two[1]]["weight"])/float(Graph[query_name][best_two[0]]["weight"])<0.6:
            best_two = [best_two[0]]
    for end  in best_two:
        New_Graph[node][end]=""
        start,end = reverse( [node,end] )
        New_Graph[start][end]=""


for line in open(END.name,'rU'):
    line = line.strip()
    line_l = line.split()
    frame =  int(line_l[8])
    query_name = line_l[0]+"+"
    subject_name= line_l[1]
    if frame==1:
        dir="-"
    else:
        dir="+"

    subject_name +=dir
    if int(line_l[5])<25:
        cache = query_name
        query_name = subject_name
        subject_name = cache
    if query_name not  in  Graph or not len(Graph[query_name]):
        continue
    if query_name in Graph and subject_name in Graph[query_name]:
        END2.write(line+'\n')
    #best = sorted( Graph[query_name],key = lambda x: Graph[query_name][x]["weight"] )[0]
    #new_name = best.translate(libary)
    #start_name  = query_name.translate(libary)
    #remove_list = []
    #for edge in Graph[new_name]:
        #if edge != start_name:
            #remove_list.append(edge)
    #for edge in remove_list:
        #remove_bi_edge(new_name, edge)
    #if subject_name == best:
        #END2.write(line+'\n')
        #remove_list = []
        #for edge in Graph[query_name]:
            #if edge!=best:
                #remove_list.append(edge)   
        #for edge in remove_list:
            #remove_bi_edge(query_name, edge)      

