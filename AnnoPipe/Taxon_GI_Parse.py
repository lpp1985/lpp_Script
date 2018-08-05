#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/6/9
"""


from lpp import *
from collections import namedtuple
from Dependcy import *
from sqlobject import *
from optparse import OptionParser
config_hash= Config_Parse()
user = config_hash["DB"]["user"]
password = config_hash["DB"]["password"]
ip = config_hash["DB"]["ip"]
port =  config_hash["DB"]["port"]
mysql_connection = "mysql -h %s -u%s -p%s --port=%s --local-infile=1 Taxon "%(ip,user,password,port)
mysql_build = "mysql -h %s -u%s -p%s  --port=%s  --local-infile=0 "%(ip,user,password,port)
connection_string = 'mysql://%s:%s@%s:%s/Taxon'%(user,password,ip,port)    
connection = connectionForURI(connection_string)
sqlhub.processConnection = connection
class Taxon_GI(SQLObject):
    class sqlmeta:
        table="GI_Taxon"
    GI = IntCol (length=50,unique=True)
    Taxon = IntCol (length=50)
    
    taxon_index= DatabaseIndex(Taxon)
    gi_index= DatabaseIndex(GI)
    
class TaxonName(SQLObject):
    class sqlmeta:
        table="TaxonName"
    Taxon = IntCol(unique=True)
    Name = StringCol()
class Taxon_Classification(SQLObject):
    class sqlmeta:
        table="Taxon_Class"
    Taxon = IntCol (length=50,unique=True)
    Class = StringCol(length=50,unique=False)
    taxon_index= DatabaseIndex(Taxon)
    class_index =  DatabaseIndex(Class)
    
def get_or_create(model, **kwargs):
    ''' use it to get or create object from a table '''

    instance = model.selectBy(**kwargs)
    
    if instance.count():
        return instance
    else:
        instance = model(**kwargs)

        return instance
if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] 
         parse eggNOG data
   
         '''
    parser = OptionParser(usage =usage )    
    
    parser.add_option("-T", "--GI", action="store",
                      dest="gi_taxon",
                      help="gi_taxon mapping infomation")
    
    
    parser.add_option("-c", "--Class", action="store",
                      dest="classfi",
                      help="Taxon classfi")
   
    parser.add_option("-n", "--Name", action="store",
                      dest="Name",
                      help="Name Dump File")   
    
    (options, args) = parser.parse_args()
    os.system(mysql_build+"-e 'create database Taxon;'")
    os.system(mysql_build+"-e 'drop database Taxon;'")
    os.system(mysql_build+"-e 'create database Taxon;'")
    Taxon_GI.createTable(ifNotExists=True)
    Taxon_Classification.createTable(ifNotExists=True)
    TaxonName.createTable(ifNotExists=True)
    NAME_CACHE = open("Taxon_Name.cache",'w')
    NAME = open(options.Name,'rU')
    for line in NAME:
        if "scientific name" in line:
            line_l = line.split("\t|\t")
            NAME_CACHE.write( line_l[0]+'\t'+line_l[1]+'\n')
    NAME_CACHE.close()
    load_des_script = """-e 'load data local infile   "%s" into table TaxonName (taxon, name);'"""%(NAME_CACHE.name)
    os.system( mysql_connection+load_des_script   )  
    os.remove( NAME_CACHE.name )
    load_des_script = """-e 'load data local infile   "%s" into table GI_Taxon (g_i, taxon);'"""%(options.gi_taxon)
    os.system( mysql_connection+load_des_script   )    
    
    load_des_script = """-e 'load data local infile   "%s" into table Taxon_Class ( taxon,class );'"""%(options.classfi)
    os.system( mysql_connection+load_des_script   )    
    
