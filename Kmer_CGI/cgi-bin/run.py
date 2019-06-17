#!/usr/bin/python
#coding:utf-8
# -*- coding: UTF-8 -*-
from pymongo import MongoClient
from  jinja2 import FileSystemLoader,Environment
import cgi,uuid,os,jinja2,sys
reload(sys)  
sys.setdefaultencoding('utf8')  
client = MongoClient('192.168.31.82', 27017)
form = cgi.FieldStorage()
org = form.getvalue('org')

db = client[org]
start_db = db.Start
delete_db = db.Delete
related_db = db.Related

begin_list = start_db.find({}).sort("attr")

deleted_nodes = delete_db.find()


#起始作图点
start_node = form.getvalue('start_node')
if  not start_node:
    start_node = start_db.find_one()["name"]
direction = start_db.find_one({"name":start_node})["attr"]
all_nodes = related_db.find({"attr":start_node})

#撤销删除的节点
recov_nodes_list = form.getvalue('Recover')
if type(recov_nodes_list)==str:
    recov_nodes_list = [recov_nodes_list]
if type(recov_nodes_list)==list:
    for e_nodes in recov_nodes_list:
        delete_db.delete_one({"name":e_nodes})


#添加作废的节点
deleted_list = form.getvalue('deleted')
if type( deleted_list ) == str:
    deleted_list = [deleted_list]
if type( deleted_list ) ==list:
    for e_node in deleted_list:
        if delete_db.find({"name":e_node}).count()<1:
            delete_db.insert( {"name":e_node} )



#开始画图

uuid_str = str(uuid.uuid4())
TEMP_LIST = open(  uuid_str+".tsv",'w')
TEMP_LIST.write(start_node+"\n")
TEMP_LIST.close()
output = "../Graph/%s.fa"%(uuid_str)
graph = "../Graph/%s.svg"%(uuid_str)
step = form.getvalue('step')
if not step:
    step ="4"
os.system(  "./Kmer_Graph_Walk -f {output} -i {org} -j {input}  -k 41 -l {input} -r {output} -s {step}".format(input= TEMP_LIST.name,output = output,step=step,org=org  ) )

#print(  "./Draw_Kmer.py  {input}  {org} {graph}".format( input=output,org=org,graph=graph ))
now_path = os.getcwd()
os.system( """  ./Draw_Kmer.py  {input}  {org} {graph}  """.format( input=output,org=org,graph=graph ) )



template_path="../Template/"
templeloader = FileSystemLoader(template_path)
env = Environment(loader = templeloader)
template = env.get_template('Run.html')





result= template.render(
    {
        "org":org,
        "begin_list":begin_list,
        "deleted_nodes":deleted_nodes,
        "all_nodes":all_nodes,
        "graph":graph,
        "fasta":output

        }
)
print(result)
