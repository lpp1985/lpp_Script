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
    
    parser.add_option("-i", "--input", action="store",
                      dest="input",
                      help="gi_taxon mapping infomation")
    parser.add_option("-o", "--output", action="store",
                      dest="output",
                      help="output")
    
    

    
    (options, args) = parser.parse_args()
    
    END = open(options.output,'w')
    RAW = open(options.input,'rU')
    RAW.next()
    for line in RAW:
        line_l = line.split("\t")
        spe = re.search("order\s+([^\t]+)", line)
	if spe:
	    spe= spe.group(1)
	else:
	    spe = re.search("species\s+([^\t]+)",line).group(1)
        gi_instance = TaxonName.selectBy(Name=spe)
	print( spe,gi_instance.count() )
        if gi_instance.count():
            taxon = gi_instance[0].Taxon
            END.write("%s\t%d\n"%( line_l[0],taxon ))
        
