#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2018/7/25
"""

import unittest
from lpp import *
from networkx import DiGraph
import string
libary=string.maketrans('+-','-+')
def get_data(RAW):
    graph=DiGraph()
    for line in RAW:
        line_l = line.strip().split("\t")
        graph.add_edge(line_l[0],line_l[1])
    return  graph
def GetPred(graph,node):
    result = ""
    if node not in graph:
        return result
    else:
        all_p = graph.predecessors(node)
        for e_p in all_p:
            return e_p
    return result
if __name__ == '__main__':
    RAW = open(sys.argv[1],'rU')
    RAW2 =  open(sys.argv[2],'rU')
    END = open("Combined_Network.tsv",'w')
    old_graph =get_data(RAW)
    ref_graph = get_data(RAW2)
    need_pred = []
    need_succ = []
    add_rela = []
    for node in old_graph.nodes():
        pred = GetPred(old_graph, node)

        
        if len(pred)==0:
            need_pred.append(node)

    
    for e_p in need_pred:
        print("Start is %s"%(e_p))
        if e_p in ref_graph:
            pred = GetPred(ref_graph,e_p)
            rev_pred = pred.translate(libary) 
            while ( pred not in old_graph and rev_pred not in old_graph) or( (pred in old_graph and  not old_graph[pred]) or ( rev_pred in old_graph and not GetPred(old_graph,rev_pred) )):
                if pred =="":
                    break
                old_graph.add_edge( pred,e_p )

                e_p = pred
                pred = GetPred(ref_graph,e_p)
                

    for node in old_graph.nodes():

        if len(old_graph[node]) ==0:
            need_succ.append(node)

    
            
    for e_s in need_succ:
        if  e_s in ref_graph:
            succ = ref_graph[ e_s ]
            if len(succ) ==0:
                break
            succ=succ.keys()[0]
            rev_succ = succ.translate(libary) 
            while (succ not in old_graph and rev_succ not in old_graph )or(       (rev_succ in old_graph  and not old_graph[rev_succ]) or(  succ in old_graph and not GetPred(old_graph,succ)  )     ):
                
                old_graph.add_edge( e_s,succ )
                #print(e_s,succ)
                e_s = succ
                succ = ref_graph[ e_s ]       
                if len(succ)==0:
                    break
                succ = succ.keys()[0]
    for start ,end in ref_graph.edges():
        if start not in old_graph and end not in old_graph and start.translate(libary) not in old_graph and end.translate(libary) not in old_graph:
            old_graph.add_edge( start,end )
 
            
    need_pred=[]
    i=0
    for node in old_graph.nodes():
        pred = GetPred(old_graph, node)

        
        if len(pred)==0:
            need_pred.append(node)   
    has = {}
    for e_p in  need_pred:
        END.write("Scaffold%s\t%s"%(i,e_p))
        has[e_p] = ""
        i+=1
        succ = old_graph[ e_p ]
        while len(succ)==1:
            e_p =  succ.keys()[0]
            if e_p in has:
                break
            END.write('; '+ e_p  )
            
            has[e_p]=""
            succ = old_graph[ e_p ]
        END.write("\n")
        
        
    print(len(has))
    