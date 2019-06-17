#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2019/4/9
"""

from pymongo import MongoClient
from  jinja2 import FileSystemLoader,Environment
import cgi,uuid,os,jinja2,sys
reload(sys)  

template_path="../Template/"
sys.setdefaultencoding('utf8')  
client = MongoClient('192.168.31.82', 27017)
form = cgi.FieldStorage()
db = {}
for i in client.database_names():
    if i in ["admin","test","local"]:
        continue
    else:
        db[i]=""
templeloader = FileSystemLoader(template_path)
env = Environment(loader = templeloader)
template = env.get_template('Index.html')
result= template.render(
        {
                "sample_list":db
            

        }
)
print(result)

