#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/6/6
"""
from  termcolor import colored
from copy import deepcopy
from optparse import OptionParser
import subprocess
import shlex
from lpp import *
from  jinja2 import FileSystemLoader,Environment
import json
import networkx as nx
import numpy as np
from networkx.readwrite import json_graph
from networkx.algorithms.clique import find_cliques
import string
libary=string.maketrans('+-','-+')
Assembled_Contig_no = 1
def reverse_path(data):
    #转换成为反向节点
    return '; '.join( [x.translate(libary) for x in data[::-1] ] )

class  Contig_Graph(nx.DiGraph):
    def __init__(self, data=None, **attr):
        super(Contig_Graph,self).__init__(data,**attr)
        self.new_id = 0
        self.assembly_detail = {}
    @staticmethod
    def reverse(data):
        #转换成为反向节点
        return [x.translate(libary) for x in data[::-1] ]

    def add_bi_edge(self,start,end):

        super(Contig_Graph,self).add_edge(start,end)
        start,end = self.reverse( [start,end] )

        super(Contig_Graph,self).add_edge(start,end)

    def add_bi_cycle(self,path):
        self.add_cycle(path)
        path_rev = self.reverse(path)
        self.add_cycle(path_rev)
    def remove_bi_edge(self,start,end):
        if self.has_edge(start,end):
            super(Contig_Graph,self).remove_edge(start,end)
        start,end = self.reverse([start,end])
        if self.has_edge(start,end): 
            super(Contig_Graph,self).remove_edge(start,end)
    def remove_bi_edges_from(self,edges_list  ):
        for key in edges_list:
            self.remove_bi_edge(key)

    def add_bi_edges_from(self,edges_list   ):
        for key in edges_list:
            self.add_bi_edge(key)

    def remove_bi_path(self,paths):
        for i  in xrange(1,len(paths)):
            if self.has_edge(paths[i-1],paths[i]):
                self.remove_bi_edge(paths[i-1],paths[i])
    def add_bi_path(self,paths):
        for i  in xrange(1,len(paths)):
            if paths[i] not in self.node:
                self.add_node(paths[i],name = paths[i])
                self.add_node(paths[i].translate(libary),name = paths[i].translate(libary))
            self.add_bi_edge(paths[i-1],paths[i])	
    def remove_bi_node( self,node  ):
        if not node.endswith('+') and not node.endswith("-"):
            node = node+'+'
        self.remove_nodes_from( [node, node.translate(libary)]   ) 
    def add_bi_node(self,node):
        if not node.endswith('+') and not node.endswith("-"):
            node = node+'+'
        self.add_nodes_from( [node, node.translate(libary)]   )  

    def unique_pred(self,start):
        if len(self.predecessors(start ))==1 :
            predes = self.predecessors(start )[0]
            if len(self.successors(predes ))==1:
                return predes
            elif len(self.successors(predes ))>1:
                return "Multi"
        elif len(self.predecessors(start ))>=1:
            return "Multi"
        else:
            return None
    def unique_succ(self,start):
        if len(self.successors(start ))==1 :
            succ = self.successors(start )[0]
            if len(self.predecessors(succ ))==1:
                return succ
            elif len(self.predecessors(succ ))>1:
                return "Multi"
            else:
                return None
        elif len(self.successors(start))>1:
            return "Multi"
        else:
            return None	


    def has_repeat_path(self,start,end,unique_contig):
        assert start[:-1] in unique_contig
        assert end[:-1] in unique_contig
        if start not in self.node or end not in self.node:
            return None
        No_Other_Graph = deepcopy(self) 

        for each_key in unique_contig:
            if each_key not in [ start[:-1],end[:-1]  ]:
                No_Other_Graph.remove_bi_node(each_key)
        if start.translate(libary) in self.node:
            No_Other_Graph.remove_node(start.translate(libary))
        if end.translate(libary) in self.node:
            No_Other_Graph.remove_node(end.translate(libary))

        if nx.algorithms.has_path(No_Other_Graph, start, end):
            end_path_list = []
            for path in nx.algorithms.all_simple_paths(No_Other_Graph,start,end):
                end_path_list.append(path)
            return end_path_list
        else:
            return None
    def check_path(self,start,end ,unique_contig ):
        candidate_path = self.has_repeat_path(start, end, unique_contig)
        if candidate_path:
            if len(candidate_path)==1:
                candidate_path = candidate_path[0]
                all_children= self.successors(start)
                for each_child in all_children:

                    self.remove_bi_edge(start,each_child)


                all_parents = self.predecessors(end)		
                for each_parent in all_parents:

                    self.remove_bi_edge(each_parent,end)				

                if len(candidate_path)>2:
                    candidate_path = candidate_path[1:-1]

                    if "; ".join(candidate_path) not in self.assembly_detail:
                        name = "New%s"%(self.new_id)
                        self.assembly_detail["; ".join(candidate_path)] = [name+'+',1]
                        self.assembly_detail['; '.join(self.reverse(candidate_path))] = [name+'-',1]
                        self.new_id+=1

                        MAPPING.write(
                            "%s\t%s\tDoubt\n"%(name,'; '.join(candidate_path)  
                                               )  
                        )
                        new_name = name+'+_1'
                        self.add_node(new_name,name =new_name   )
                        rev_name = new_name.translate(libary)
                        self.add_node(rev_name, name =rev_name )
                    else:
                        self.assembly_detail["; ".join(candidate_path)][1]+=1
                        self.assembly_detail['; '.join(self.reverse(candidate_path))][1]+=1
                        new_name = self.assembly_detail['; '.join(candidate_path)][0] + "_"+str(self.assembly_detail['; '.join(candidate_path)][1] )
                        rev_name = new_name.translate(libary)
                        self.add_node(new_name,name =  new_name  )
                        self.add_node(rev_name,name =  rev_name  )

                    self.add_bi_path([start,new_name,end])
                else:
                    self.add_bi_edge(start, end)
            else:
                print(candidate_path)

    def same(self,raw,new):
        for each_succ in self.successors(raw):
            self.add_bi_edge(new,each_succ)
        for each_pre in self.predecessors(raw):
            self.add_bi_edge(each_pre, new)
def Remove_rev_Singleton(  G   ):
    '''将所有的反向singleton删除，这样得到的图比较好看'''
    singleton_deleted = []
    for node in G:
        if  not G.successors(node) and not G.predecessors(  node  ):
            if node.endswith('-'):
                singleton_deleted.append( node )
    G.remove_nodes_from( singleton_deleted  )	
def Alle_Removal(G):
    deleted = {}
    for eachNode in G.nodes():
        succers = G.successors(eachNode)
        if len(succers)==2:
            eachsucc1,eachsucc2 = succers
            if eachsucc1 not in deleted and eachsucc2 not in deleted: 
                if len(  G.predecessors(eachsucc1) )==1 and len(  G.predecessors(eachsucc2) )==1 and G.successors(eachsucc1) == G.successors(eachsucc2)  and len(  G.successors(eachsucc1)  )==1:
                    if eachsucc1[:-1] != eachsucc2[:-1]:
                        deleted[eachsucc1] = ""
                        deleted[ eachsucc1.translate(libary) ] = ""
    New_Graph = deepcopy(G)

    for node in deleted:
        New_Graph.remove_bi_node(node)
    return New_Graph,deleted


def Orphan_Removal( G ):
    '''删除一系列寡妇节点，如
    1--------->2------->3------>4
                ------->5
    则删除5
    '''
    deleted = {}
    New_Graph = deepcopy(G)

    for eachNode in G.nodes():
        if eachNode in deleted:
            continue
        if len( G.successors(eachNode) )==0 or len( G.predecessors(eachNode) )==0 :

            deleted[eachNode]=""
            new_name = eachNode.translate(libary)
            deleted[new_name]=""   
    while len(deleted)>0:
        for key in deleted:

            if key in New_Graph.node:
                New_Graph.remove_bi_node(key)
        deleted = {}
        G = New_Graph
        New_Graph = deepcopy(G)
        for eachNode in G.nodes():
            if eachNode in deleted:
                continue
            if len( G.successors(eachNode) )==0 or len( G.predecessors(eachNode) )==0 :
        
                deleted[eachNode]=""
                new_name = eachNode.translate(libary)
                deleted[new_name]=""         
            
    return New_Graph

def Transitive_Remove( G ):
    '''对于一些tiling结构进行修剪，遍历下游节点，如果下游节点的下游节点与上上游节点相连，
    则将该边删除
    '''
    G_output = deepcopy(G)
    isinstance( G_output,Contig_Graph  )
    isinstance(    G,Contig_Graph     )
    for node in G:
        if len( G.successors( node ) )==2:

            for each_child in G.successors(  node ):
                for each_child_child in G.successors_iter( each_child  ):
                    if G.has_edge( node,each_child_child    ):
                        G_output.remove_bi_edge(  node,each_child_child    )

    return G_output
def Get_Long_edges( G ):
    '''对于一些tiling结构进行修剪，遍历下游节点，如果下游节点的下游节点与上上游节点相连，
    则将该边删除
    '''
    G_output = deepcopy(G)
    isinstance( G_output,Contig_Graph  )
    isinstance(    G,Contig_Graph     )
    for node in G:
        if len( G.successors( node ) )>1:

            for each_child in G.successors(  node ):
                for each_child_child in G.successors_iter( each_child  ):
                    if G.has_edge( node,each_child_child    ):
                        G_output.remove_bi_edge(  node,each_child_child    )

    return G_output

def Get_list(  has_list,location_list  ):
    #将对应的location的元素打入到输出列表中
    output_list = []
    for i in location_list:
        output_list.append( has_list[ i ]   )
    return output_list

def Get_Contig( G ,DETAIL,name_prefix ):
    global Assembled_Contig_no
    Assembled_Contig_no = 1
    """获得contig，并且将一些单链和成环的contig删除"""
    #G_Trimed输出被检查过的结果
    #has_checked输出不需要检查的节点（已经是其他contig的子路径，或者该路径已经被查询完毕）
    #deleted_node输出冗余和检查过的双向图结果，即正向已经被证明是singleton或者plasmid，反向直接扔掉
    #另外，deleted_node还记录已经被检查完，放入contig的reads，在最后的图中，统一删除
    #Contig_relation 记录每一个contig的连接关系

    all_path = {}
    deleted_node = {}
    plasmid_contig = {}
    Contig_relation = {}
    G_Output =  Contig_Graph()
    rela_5 = {}
    rela_3 = {}
    #DETAIL记录Contig的明细信息



    #这是个计数器，用来最后contig1的编号 

    isinstance(G,Contig_Graph)
    def __Enlongation_5(node):
        in_path = {node:""}
        node = G.unique_pred(node)
        status = node
        while node !="Multi":
            if node :
                if node not in in_path:
                    path.insert(0,node )
                    in_path[node] = ''
                    node = G.unique_pred(node)
                    status = node
                else:
                    status ="Circle"

                    break
            else:
                break

        return status

    def __Enlongation_3(node):
        node = G.unique_succ(node)
        status_3 = node
        while node !="Multi":
            if node :
                path.append(node )
                node = G.unique_succ(node)
                status_3 = node

            else:
                break

        return status_3	


    #def __Enlongation_5(node):
        ##查询overlap图，向5’端延伸，删除成环的节点
        #if node in has_checked :
            #return 1
        #already_in_path[ node ] = ''
        #has_checked[ node ] = ''

        #if len(  G.predecessors(node)  ) ==1:
            #predecessor  = G.predecessors(node)[0]
            #if  len( G.successors( predecessor ) )==1:
                #if predecessor not in already_in_path:
                    #already_in_path[ predecessor  ] = ''
                    #path.insert( 0 ,predecessor)

                    #status = __Enlongation_5( predecessor)
                    #return status
                #else:
                    #status =  2
                    #return status

            #else:
                #status = 1
                #return status

        #elif  len(  G.predecessors(node)  )==0:
            #status =0
            #return status
        #else:
            #status =1
            #return status
    #def __Enlongation_3(node):
        ##查询overlap图，删除成环的节点
        #global status
        #has_checked[node] = ''

        #if len(G.successors(node))<2:

            #successor  = G.successors(node)

            #if successor:
                #successor = successor[0]

                #if successor in has_checked:
                    #status = 1
                    #return status

                #if  len( G.predecessors( successor ) )==1:

                    #already_in_path[ successor  ] = ''
                    #path.append( successor)

                    #status = __Enlongation_3( successor)
                    #return status

                #else:
                    #status = 1
                    #return status
        #else:
            #status =1	
            #return status


    def __Check_result():


        """检查最后的path结果"""
        global Assembled_Contig_no
        new_contig_name =  "%s%s" % ( name_prefix,  Assembled_Contig_no )

        Assembled_Contig_no+= 1
        #迁移连接关系
        #如果status==1，则该contig两边都连接不确定的关系，该contig的正反向都保留
        #如果status==2，则该contig直接就是一个质粒，所有的正向、反向节点都删除，
        #如果status==0，则该contig为一个链状contig，删除反向节点，不做关系迁移
        #所有的相关节点都准备删掉
        Contig_relation[new_contig_name] = path
        for d_node in path:
            traversed_node =  d_node.translate(libary)
            deleted_node[ traversed_node ]  = ''
            deleted_node[ d_node ]  = ''
        ##链状和环状contig去除
        if status !="Circle":
            rela_5[path[0]]=new_contig_name+'+'
            rela_3[new_contig_name+'+']=path[-1]
            if  path[-1].translate(libary) in G.node:
                rela_5[ path[-1].translate(libary)]=new_contig_name+'-' 
            if path[0].translate(libary)  in G.node:
                rela_3[new_contig_name+'-']=path[0].translate(libary)			
        if  status:
            #删除节点


            if status == "Multi":
                G_Output.add_node(new_contig_name+'+', name = new_contig_name+'+')
                G_Output.add_node(new_contig_name+'-', name = new_contig_name+'-')
                category = "Doubt"
            else:
                category =  "Plasmid"
                plasmid_contig[new_contig_name]  =  path
                G_Output.add_edge(new_contig_name+'+',  new_contig_name+'+') 
            DETAIL.write (new_contig_name+'\t'+'; '.join(path)+"\t"+category+'\n'  )

        else:

            category =  "Chain"
            #节点迁移
            
            DETAIL.write (new_contig_name+'\t'+'; '.join(path)+"\t"+category+'\n'  )
            G_Output.add_node(new_contig_name+'+', name = new_contig_name+'+' )
            G_Output.add_node(new_contig_name+'-', name = new_contig_name+'-' )
        all_path[new_contig_name+'+'] = path






    for each_node in G:
        if each_node  in deleted_node:
            continue
        #path用来计路径，already_in_path用来记录那些节点已经出现在路径中
        #status用来记录是否需要将该路径对应的反向节点删除
        path = [each_node]


        #status 可以为三种情况，0的话是singlton节点，1的时候是开叉分支，2的时候是plasmid节点
        already_in_path = {}
        status = __Enlongation_5(each_node)



        if status!='Circle':
            status3 = __Enlongation_3(each_node)

            if status:
                status = status3


        __Check_result()
    for each_node in G_Output.node:
        if each_node in rela_3:
            for each_key in G.successors( rela_3[each_node] ):
                if each_key in rela_5:
                    G_Output.add_bi_edge(each_node,rela_5[each_key])

    return  Contig_relation, plasmid_contig,G_Output,all_path


def Get_Reference_Graph(  FILE  ):
    #读参考图文件并生成一个双向overlap图，该图记录参考序列
    Reference_graph = Contig_Graph()

    name = {}
    all_chr= []
    all_ref = []
    for line in FILE:

        line_l =line.strip().split('\t')
        name_2,reads = line_l
        if name_2 not in name:
            if len(  name ):
                all_chr.append( all_ref  )
                name[ name_2 ] = ''
                all_ref = []
        #判断上下游是否在reference上，都不在的话，跳过！
        all_ref.append( reads  )
    else:
        all_chr.append( all_ref  )        
    for each_chr in all_chr:	
        if not len(each_chr):
            continue
        Reference_graph.add_bi_cycle(each_chr)
    return Reference_graph



def Get_Overlap_Graph(FILE):
    '用best.edges进行overlap graph的构建'
    '并且记录所有的best reads的ID'
    Overlap_Graph = Contig_Graph()
    all_nodes = {}
    for line in FILE:
        if '#' in line:
            continue
        data = []
        line_l = line[:-1].split('\t')
        all_nodes[ line_l[0]  ] = ''
        all_nodes[ line_l[2] ] = ''
        all_nodes[ line_l[4]  ] = ''
        if line_l[3] =="5'":
            tag = '-'
        else:
            tag = '+'

        data.append(  line_l[2]+tag  )
        data.append( line_l[0]+'+'  )
        if line_l[5] =="5'":
            tag = '+'
        else:
            tag = '-'
        data.append( line_l[4]+tag )
        Overlap_Graph.add_bi_path(data)
    Overlap_Graph = Transitive_Remove(Overlap_Graph)
    Overlap_Graph.remove_bi_node("0+")
    return Overlap_Graph, all_nodes



def Draw_Web(file_name,G):
    templeloader = FileSystemLoader(cele_path+"/Template")
    env = Environment(loader = templeloader)
    template = env.get_template('template.html')
    root_path = os.path.split(file_name)[0]
    if not os.path.exists(  root_path+'/css' ):
        shutil.copytree( cele_path+'/Template/css', root_path+'/css'  )

    if not os.path.exists(  root_path+'/js' ):
        shutil.copytree( cele_path+'/Template/js', root_path+'/js'  )     
    END = open( file_name ,'w' )
    web_data = json_graph.node_link_data(G)

    END.write(
        template.render(

            raw_links = json.dumps(web_data['links']    )    ,
            raw_nodes = json.dumps(web_data['nodes']   )  
        )
    )	

def Find_Overhang(Raw_queryID,Raw_subjectID,overlstore):
    """寻找overlap关系中存在的overhang信息"""
    overlap_PATH =celera_assembler+'/overlapStore'

    if Raw_queryID.endswith("+"):
        direction = '3'
        queryID = Raw_queryID
        subjectID = Raw_subjectID	
    elif  Raw_subjectID.endswith('+'):
        direction = '5'
        queryID = Raw_queryID
        subjectID = Raw_subjectID
    else:
        queryID = Raw_subjectID[:-1]+'+'
        subjectID = Raw_queryID[:-1]+'+'
        direction='5'

    #检查是否是3‘或者5’端overlap

    if direction =='3':
        command =  overlap_PATH+''' -d %(overlap)s     -b %(query)s -e %(query)s -d3 | grep -P "\s+%(subject)s\s+"'''%( 
            {
                'overlap':overlstore,'query':queryID[:-1],'subject':subjectID[:-1] }  
        )


        detailout = os.popen(
            command 
        )


    else:
        command =  overlap_PATH+''' -d %(overlap)s     -b %(subject)s -e %(subject)s -d5 | grep -P "\s+%(query)s\s+"'''%( 
            {
                'overlap':overlstore,'query':queryID[:-1],'subject':subjectID[:-1] 
            }  
        )



        detailout = os.popen(
            command 
        )
    detail_block = None
    for line in detailout:
        detail_block = line.strip().split()
        break
    #detail_block = detailout.next().strip().split()

    if not detail_block:

        print("error!!!")
        return ""

    else:
        if Raw_queryID[-1]=="+" and Raw_subjectID[-1]=="-":
            location = -int(detail_block[4])
        elif Raw_queryID[-1]=="-" and Raw_subjectID[-1]=="-":
            location = int(detail_block[3])
            if len(str(location))>6:
                location = location - int("100000001",16)
        elif Raw_queryID[-1]=="+" and Raw_subjectID[-1]=="+":
            location = -int(detail_block[4])
        elif  Raw_queryID[-1]=="-" and Raw_subjectID[-1]=="+":
            location = int(detail_block[4])
            if len(str(location))>6:
                location = location - int("100000001",16)			
    return location


def Mummer_parse(  file_name  ):
    """运行nucmer并且根据show-coords的数据对所有的结果进行解析
    得到Contain，BEGIN，END，Identity的数据并将自身比自身的序列去除。而后，将这些数据分门别类放到一个多维哈希
    output_category中，其中一维键为tag，二维键是行

    """
    #运行nucmer
    output_category = Ddict()
    cache_path = os.path.dirname(os.path.abspath(file_name))
    os.system(  "nucmer --maxmatch   %(genome)s   %(genome)s   -p %(path)s/cache >/dev/null 2>&1"%( {"genome":file_name,"path":cache_path }  )  )
    contain_data =  os.popen(  """ delta-filter  -i 95 %s/cache.delta > %s/filter.delta&&show-coords  %s/filter.delta -odTlb -L 40  | grep -P  "\[\CONTAINED\]$" """%( cache_path, cache_path,cache_path )      )
    already_contained = {}
    for line in contain_data:
        line_l = line.split("\t")
        contained = line_l[-3]
        already_contained[contained] = ''
    align_data = os.popen(  """show-coords  %s/cache.delta -odTl -L 40  | grep -P  "\[\S+\]$" """%( cache_path )      )
    align_data.next()
    already_contained = {}
    for line in align_data:
        line_l = line.strip().split("\t")
        identy,refname,queryname,tag  = Get_list( line_l,[6,11,12,13]   )
        if refname ==queryname:
            continue
        if queryname in already_contained or refname in already_contained:
            continue
        if float(  identy  )<90:
            continue
        name = re.search( "\[(\w+)\]",tag  ).group(1)
        output_category[ name  ][ line ] = ''

    return output_category


def Relation_parse(  file_name ,output_category,name_prefix,contain_trim=1   ):
    """
    接受所有的ContigID
    根据Mummer_parse获得的结果，对identity的节点，和contain的节点进行删除，取长的保留。
    而后将保留的节点根据Begin和END的关系构建overlap图
    对获得的Contig图进行修订用_Trim函数
    该函数的作用如下：
    如果 a-->b b-->c 并且 a-->c
+          ----------------------------------
    -------------------
             -----------------------------------------------
    这种情况会导致双向唯一过敏感，
    纠正的方法是：
    如果某一个节点的下游开叉，
    遍历开叉下游节点，对开叉下游节点的连接关系进行遍历，如果下下游节点存在于上上游节点的下游关系中，则将该关系删除






    """













    def _Trun_contain_Rela(G_start,contain_rela):

        G = deepcopy(G_start)
        for big in contain_rela:
            if big in  G_start.succ:
                big_succ = set(G_start.succ[big])
            else:
                big_succ=set([])

            if big in G_start.pred:
                big_pred =set( G_start.pred[big] )
            else:
                big_pred =set([])
            for small in contain_rela[big]:

                if small in  G_start.succ:
                    small_succ = set(G_start.succ[small] )
                else:
                    small_succ=set([])

                if small in G_start.pred:
                    small_pred = set( G_start.pred[small] )
                else:
                    small_pred =set([])	
                #print("Big is %s!!"%(big))
                #print("Small is %s!!"%(small))
                #print("Big succ is %s!!"%(big_succ))
                #print("Small succ is %s!!"%(small_succ))
                #print("big pred is %s!!"%(big_pred))
                #print ("Small pred is %s!!"%(small_pred))
                    
                if small_succ<=big_succ and small_pred <=big_pred:
                    #print("Remove %s"%(small))
                    G.remove_bi_node(small)
                   
            #print("###############")
                #if big =="XV05_100+":
                    #print(big,small)
                    #print(big_succ,small_succ)
                    #print(big_pred,small_pred)


        return G



    #检查contain关系，并将所有contain的节点删除

    containeed_line = {}

    #首先对所有的reads做筛选，将contain和完全相同的去除掉

    #检查contain关系，并将所有contain的节点删除
    CONTAIN = open( name_prefix+'.contains','w'  )
    containeed_line = {}
    CONTAIN.write( "Contig\tIncluded_by\n"  )
    #从图中去除所有的contain关系，如果这两个关系的上下游都一样，则讲这个节点替换为大节点,contain_rela记录所有的包含关系，这类关系最后通过图
    #的结构进行削减，如果小的overlap关系和大的overlap关系完全相同，则将该节点删除


    contain_rela = Ddict()
    for each_line in output_category[ "CONTAINS"  ]:
        line_l = each_line.strip().split("\t")
	
        big_contig,small_contig = Get_list( line_l,[11,12] )
        frame = line_l[10]
        if frame =="-1":
            frame ="-"
        else:
            frame ='+'
        contain_rela[big_contig+'+'][small_contig+frame] = ""
        CONTAIN.write( small_contig+'\t'+big_contig+'\n'   )
    #检查contained关系
    for each_line in output_category[ "CONTAINED"  ]:
        line_l = each_line.strip().split("\t")
        big_contig,small_contig = Get_list( line_l,[12,11] )
        frame = line_l[10]
        if frame =="-1":
            frame ="-"
        else:
            frame ='+'		
        contain_rela[big_contig+frame][small_contig+"+"] = ""
        CONTAIN.write( small_contig+'\t'+big_contig+'\n'   )

    for data in containeed_line:
        CONTAIN.write(data)
    #检查identity
    #print(len( output_category[ "IDENTITY"  ]))
    for each_line in output_category[ "IDENTITY"  ]:

        line_l = each_line.strip().split("\t")
        query_length,subject_length,frame,query_contig,subject_contig = Get_list( line_l,[7,8,10,11,12] )
        query_length =int( query_length );subject_length = int( subject_length )
        query_contig = query_contig+"+"
        if frame=="-1":
            frame ="-"
        else:
            frame ="+"
        subject_contig = subject_contig+frame
        if query_length>subject_length:
            small_contig = subject_contig
            big_contig = query_contig
        elif query_length==subject_length:
            big_contig,small_contig = sorted([query_contig,subject_contig]) 
        else:
            small_contig = query_contig
            big_contig = subject_contig

        contain_rela[big_contig][small_contig] = ""
        CONTAIN.write( small_contig+'\t'+big_contig+'\n'   )			
    #print(contain_rela) 
    #开始对 Start和END进行解析，期间只要碰到 删除的节点，就跳过
    #设置一个初始图
    G_contig_init = Contig_Graph()
    # 添加所有的节点 添加双向
    FASTA = fasta_check(open( file_name ,'rU') )
    for t,s in FASTA:
        each_contig = t.strip()[1:].split()[0]

        for direction in ['-','+']:
            G_contig_init.add_node( each_contig+direction,name = each_contig+direction )

    #根据Begin和END添加关系
    #根据BEGIN进行构图
    for each_line in output_category[ "BEGIN"  ]:
        line_l = each_line.split("\t")
        query_node,subj_node,frame = Get_list( line_l,[ 11,12,10 ]  )
        #如果删除的节点出现了这些数据则跳过
        #否则
        if frame =="1":
            direction = "+"
        else:
            direction = '-'

        G_contig_init.add_bi_edge(  subj_node+direction, query_node+'+'    )




    #根据END进行构图
    for each_line in output_category[ "END"  ]:
        line_l = each_line.split("\t")
        query_node,subj_node,frame = Get_list( line_l,[ 11,12,10 ]  )
        #如果删除的节点出现了这些数据则跳过
        #否则
        if frame =="1":
            direction = "+"
        else:
            direction = '-'

        G_contig_init.add_bi_edge( query_node+'+', subj_node+direction    )







    #修图，将tiling边删除
    #print(contain_rela)
    if contain_trim==1:
        G_contig_init = _Trun_contain_Rela(G_contig_init, contain_rela)
 
    G_Trimed = Transitive_Remove( G_contig_init )
    #G_Trimed = G_contig_init


    #双向图变成单向图
    ##首先，将singleton的反向节点删除
    Remove_rev_Singleton(G_Trimed)

    #得到新的最终的Graph


    return G_Trimed
    ##到目前为止 G_Trimed完成，下面遍历所有的contig，遍历标准是predecessor 为1，successors为1，认为是一个contig
    ##见这些contig的明细信息保存下来



def Find_Cycle(Graph):
    #寻找所有可能的边
    all_cycle = []
    already_cycle = []
    for each_cycle in nx.simple_cycles(Graph):

        this_graph =nx.DiGraph()
        this_graph.add_cycle(each_cycle)
        edges = this_graph.edge
        if edges not in already_cycle:

            already_cycle.append(edges)
            cycle_rev = [x.translate(libary) for x in each_cycle[::-1] ]
            this_graph =nx.DiGraph()
            this_graph.add_cycle(cycle_rev)
            already_cycle.append( this_graph.edge)

            all_cycle.append(each_cycle)
        else:
            continue

    return all_cycle

def Reverse_Path(path):
    rev_path = []
    for node in path[::-1]:
        rev_path.append( node.translate(libary) )
    return rev_path









