#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/6/20
"""
from Dependcy import *
from ConfigParser import ConfigParser
import os,sys,redis
sys.path.append(
    os.path.abspath(os.path.split(__file__)[0]+'../Lib/' )
)	
from lpp import *
from collections import namedtuple

from sqlobject import *
from optparse import OptionParser
class NOG_des(SQLObject):
    class sqlmeta:
        table="NOG_Description"
    Name = StringCol(length=50,unique=True)
    Description = StringCol()
    name_index= DatabaseIndex(Name)

class NOG_GENE( SQLObject  ):
    class sqlmeta:
        table="Gene_NOG"
    Gene = StringCol(length=100)
    NOG = StringCol(length=50,)
    gene_index = DatabaseIndex(Gene)
    nog_index = DatabaseIndex(NOG)
    
    
class NOG_CAT( SQLObject  ):
    class sqlmeta:
        table="NOG_Category"
    NOG =StringCol(length=50)
    Cat = StringCol(length=10)
    nog_index = DatabaseIndex(NOG)
    cat_index = DatabaseIndex(Cat)
    
class CAT_DES( SQLObject  ):
    class sqlmeta:
        table="CAT_des"
    Abb= StringCol(length=4)
    Description = StringCol()
    cat_index = DatabaseIndex(Abb)
    
        




def get_or_create(model, **kwargs):
    ''' use it to get or create object from a table '''

    instance = model.selectBy(**kwargs)
    
    if instance.count():
        return instance
    else:
        instance = model(**kwargs)

        return instance
config_hash= Config_Parse()
user = config_hash["DB"]["user"]
password = config_hash["DB"]["password"]
port =  config_hash["DB"]["port"]
ip = config_hash["DB"]["ip"]
mysql_connection = "mysql -h %s -u%s -p%s --port=%s  --local-infile=1 eggNOG "%(ip,user,password,port)
mysql_build = "mysql -h %s -u%s -p%s --port=%s --local-infile=1 "%( ip, user, password, port )
connection_string = 'mysql://%s:%s@%s:%s/eggNOG'%( user,password,ip,port )    
connection = connectionForURI(connection_string)
sqlhub.processConnection = connection

if __name__ == '__main__':
      
    
    usage = '''usage: python2.7 %prog [options] 
         parse eggNOG data
   
         '''
    parser = OptionParser(usage =usage )    
    parser.add_option("-d", "--DES", action="store",
                      dest="description",
                      help="NOG descripton")
    parser.add_option("-m", "--Mapping", action="store",
                      dest="mapping",
                      help="Gene NOG Mapping")
    parser.add_option("-a", "--Anno", action="store",
                      dest="anno",
                      help="Cat annotation")    
    parser.add_option("-c", "--Cat", action="store",
                      dest="cate",
                      help="NOG gene Category")  


    (options, args) = parser.parse_args()
    
    os.system(mysql_build+"-e 'create database eggNOG;'")
    os.system(mysql_build+"-e 'drop database eggNOG;'")
    os.system(mysql_build+"-e 'create database eggNOG;'")
    general_config = ConfigParser()
    path = os.path.split(os.path.abspath(__file__))[0]+'/'
    general_config.read(
        os.path.join( path+"database.ini")
    ) 


    NOG_des.createTable(ifNotExists=True)
    NOG_GENE.createTable(ifNotExists=True)
    CAT_DES.createTable(ifNotExists=True)
    NOG_CAT.createTable(ifNotExists=True)

    ANNO = open(options.anno,'rU')


    for (abb,description) in  re.findall( "\[(\w)\]\s+([^\n]+)",ANNO.read()  ):
        get_or_create(CAT_DES,Abb=abb,Description=description)
    
    
    
   
  
  

    load_des_script = """-e 'load data local infile   "%s" into table NOG_Description (name, description);'"""%(options.description)
    os.system( mysql_connection+load_des_script   )


    load_cate_script = """-e 'load data local infile   "%s" into table  NOG_Category (no_g,cat);'"""%(options.cate)
    os.system( mysql_connection+load_cate_script   )

    load_rela_script = """-e 'load data local infile   "%s" into table Gene_NOG (gene,no_g);'"""%(options.mapping)
    os.system( mysql_connection+load_rela_script   )


   
    
    
