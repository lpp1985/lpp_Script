#!/usr/bin/python2
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/26
from lpp import *
from sqlobject import *
from optparse import OptionParser
from Dependcy import *
config_hash= Config_Parse()
user = config_hash["DB"]["user"]
password = config_hash["DB"]["password"]
port =  config_hash["DB"]["port"]
ip = config_hash["DB"]["ip"]
mysql_connection = "mysql -h %s -u%s -p%s --port=%s --local-infile=1 GO "%(ip,user,password,port)
mysql_build = "mysql -h %s -u%s -p%s --port=%s --local-infile=1 "%( ip,  user,password,port )
connection_string = 'mysql://%s:%s@%s:%s/GO'%(user,password,ip,port)    
print(connection_string)
connection = connectionForURI(connection_string)
sqlhub.processConnection = connection

class GO_ROOT(SQLObject):
    class sqlmeta:
        table="GO_ROOT"
    Go = StringCol (length=50)
    Rank = IntCol (length=2)

    go_index= DatabaseIndex(Go)
    rank_index= DatabaseIndex(Rank)

class GO_ALTER(SQLObject):
    class sqlmeta:
        table="GO_ALTER"	
    Go_raw = StringCol (length=50)
    Change_to = StringCol (length=50)
    raw_index= DatabaseIndex(Go_raw)
    change_index= DatabaseIndex(Change_to)


class GO_DEF(SQLObject):
    class sqlmeta:
        table="GO_DEF"	
    Go = StringCol (length=50,unique = True)
    Def = StringCol()
    go_index= DatabaseIndex(Go)
    
class UNIPROT_GO(SQLObject):
    class sqlmeta:
        table="UNIPROT_GO"    
    UniID = IntCol()
    Go = StringCol (length=50)
    go_index= DatabaseIndex(Go)
    uniprot_index = DatabaseIndex(UniID)
class UNIPROT(SQLObject):
    class sqlmeta:
        table="UNIPROT" 
    
    Uniprot = StringCol (length=50)
    UniID = IntCol()
    uniprot_index = DatabaseIndex(Uniprot)
    
class UNIPROT_GI(SQLObject):
    class sqlmeta:
        table="UNIPROT_GI" 
    
    
    UniID = IntCol()
    GI = IntCol()
    uniprot_index = DatabaseIndex(UniID)
    gi_index = DatabaseIndex(  GI  )
    
class GO_SON(SQLObject):
    class sqlmeta:
        table="GO_SON"
    Father = StringCol (length=50)
    Son = StringCol (length=50)
    Rela = StringCol (length=50)
    father_index= DatabaseIndex(Father)
    son_index  = DatabaseIndex(Son)
    rela_index = DatabaseIndex(Rela)
class GO_COMPONENT(SQLObject):
    class sqlmeta:
        table="GO_COMP"
    Go = StringCol (length=50)
    Compent = StringCol()
    go_index= DatabaseIndex(Go)
    
if __name__ == '__main__':
    usage = '''usage: python2.7 %prog [options] 
         parse eggNOG data

         '''
    parser = OptionParser(usage =usage )    
    parser.add_option("-I", "--OBO", action="store",
                      dest="Input",
                      help="obo file")
   
    parser.add_option("-M", "--Mapping", action="store",
                      dest="Mapping",
                      help="mapping file")     

    (options, args) = parser.parse_args()


    RAW = block_reading(  open( options.Input,'rU'   ) ,tag='^\[\S+\]'       )
    FATHER = open( 'relationship.son','w'  )
    NAME_DEF = open( 'NAME_DEF.list','w' )

    root = {}
    ALTER = open( 'relationship.alter','w'   )
    leaf_2 = {}
    leaf_3 = {}

    COMPONENT = open("GO_Component",'w')
    #os.system(mysql_build+"-e 'create database GO;'")
    os.system(mysql_build+"-e 'drop database GO;'")
    os.system(mysql_build+"-e 'create database GO;'")    
    for e_b in RAW:
        if  e_b .startswith('id: GO:'):
            go_id = re.search( '^id\: (GO\:\d+)',e_b ).group(1)
            alter_all = re.findall( 'alt_id\: (\S+)',e_b )
            for each_altered in alter_all:
                ALTER.write( each_altered+'\t'+go_id+'\n' )
            consider_all = re.findall( 'consider\: (\S+)',e_b )
            for each_consider in consider_all:
                ALTER.write( go_id+'\t'+each_consider+'\n' )


            if 'is_obsolete: true' in e_b:
                replace_by = re.findall("replaced_by\: (\S+)",e_b)
                for replace in replace_by:
                    ALTER.write( go_id+'\t'+replace+'\n' )
                continue
            compo = re.search("namespace\:\s+(\S+)",e_b).group(1)
            COMPONENT.write( go_id+'\t'+compo+'\n'   )            
            go_id = re.search( '^id\: (GO\:\d+)',e_b ).group(1)
            name = re.search( '\nname\: (.+)\n',e_b  ).group(1)
            define = re.search(  'def\: (.+)',e_b )
            if not define:
                define = ""
            else:
                define = define.group(1)
            NAME_DEF.write( go_id+'\t'+name+'\t'+define+'\n'    )
            all_father = re.findall( 'is_a\: (GO\:\d+)',e_b )
            for each_father in all_father:
                FATHER.write( each_father+'\t'+go_id+'\tis_a\n' )
            if not all_father:
                root[ go_id ] = ''
            part_of = re.findall('relationship\: part_of (\S+)',e_b)
            for each_comp in part_of:
                FATHER.write(each_comp+'\t'+go_id+'\tis_part_of\n'   )

    all_sonR = File_Ddict(  open( FATHER.name,'rU'  )  ).read(1,2)
    for each_key1 in root:
        for each_key2 in all_sonR[ each_key1  ]:
            leaf_2[ each_key2 ] = ''
            for each_key3 in all_sonR[ each_key2  ]:
                leaf_3 [ each_key3 ] = ''
    ROOT = open( 'ROOT.root','w'  )
    for key1 in root:
        ROOT.write( key1+'\t1\n'  )
    for key2 in leaf_2:
        ROOT.write( key2+'\t2\n'  )
      
 

    GO_ROOT.createTable(ifNotExists=True)
    GO_ALTER.createTable(ifNotExists=True)
    GO_DEF.createTable(ifNotExists=True)
    UNIPROT.createTable(ifNotExists=True)
    UNIPROT_GI.createTable(ifNotExists=True)
    UNIPROT_GO.createTable(ifNotExists=True)
    GO_SON.createTable(ifNotExists=True)
    GO_COMPONENT.createTable(ifNotExists=True)
    
    load_des_script = """-e 'load data local infile   "%s" into table GO_ROOT (go, rank);'"""%(ROOT.name)
    os.system(mysql_connection+load_des_script)
    
    load_des_script = """-e 'load data local infile   "%s" into table GO_DEF (go, def);'"""%(NAME_DEF.name)
    os.system(mysql_connection+load_des_script)
    
    load_des_script = """-e 'load data local infile   "%s" into table GO_ALTER (go_raw, change_to);'"""%(ALTER.name)
    
    os.system(mysql_connection+load_des_script)
    
    
    
    load_des_script = """-e 'load data local infile   "%s" into table GO_SON (father, son,rela);'"""%(FATHER.name)
    
    os.system(mysql_connection+load_des_script)
    
    
    load_des_script = """-e 'load data local infile   "%s" into table GO_COMP (go, compent);'"""%(COMPONENT.name)
    os.system(mysql_connection+load_des_script)
    os.system(  """GO_Uniprot_GI_Parse  -i %s  -g Uniprot_GO.db -u Uniprot.db -k Uniprot_GI.db"""%( options.Mapping )   )
   
    
    

    
        
    
    
    load_des_script = """-e 'load data local infile   "Uniprot_GO.db" into table UNIPROT_GO (uni_id, go);'"""
    os.system(mysql_connection+load_des_script)
    
    
    load_des_script = """-e 'load data local infile   "Uniprot_GI.db" into table UNIPROT_GI (uni_id, g_i);'"""
    os.system(mysql_connection+load_des_script)
    
    
 
    
    load_des_script = """-e 'load data local infile   "Uniprot.db" into table UNIPROT (uniprot, uni_id);'"""
    os.system(mysql_connection+load_des_script)

    os.remove("Uniprot_GI.db")
    os.remove("Uniprot.db")
    os.remove("Uniprot_GO.db")
