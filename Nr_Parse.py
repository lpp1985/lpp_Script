#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/3/16
"""

from os.path import abspath
import redis
from optparse import OptionParser
from lpp import *
from ConfigParser import ConfigParser
def Redis_trans(data_hash):
    out_data = ""
    if type(data_hash)==type(Ddict()):

        for key1 in data_hash:

            for key2,value in data_hash[key1].items():
                out_data+="""*4\r\n$4\r\nhset\r\n"""
                out_data+="$%s\r\n"%( len(key1) )
                out_data+="%s\r\n"%( key1 )
                out_data+="$%s\r\n"%( len(key2) )
                out_data+="%s\r\n"%( key2 )
                out_data+="$%s\r\n"%( len(value) )

                out_data+="%s\r\n\r\n"%( value )


    return out_data
usage = "python2.7 %prog [options]"
parser = OptionParser(usage =usage )
parser.add_option("-d", "--Database", action="store",
                  dest="DB_FILE",

                  help="Database File")


parser.add_option("-n", "--Nr", action="store",
                  dest="NR",

                  help="Nr fasta sequence")
general_config = ConfigParser()
path = os.path.split(os.path.abspath(__file__))[0]+'/'
general_config.read(
    os.path.join( path+"database.ini")
) 
db_number = general_config.get("Redis", "nr")    

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    
    path = os.path.split(os.path.abspath(__file__))[0]+'/'
    
    r = redis.Redis(host='192.168.31.71',port=6379,db=int(db_number))
    r.flushdb()
    DB_FILE = open( os.path.abspath(options.DB_FILE),'w')

    NR_ANNO_DETAIL = fasta_check(open(   os.path.abspath(  options.NR ),'rU'   ) )

    data_hash = Ddict()
    for t,s in NR_ANNO_DETAIL:
        t= t[1:-1]
        name = t.split()[0]
        data_hash[name]["Annotation"] = t
        s1 = re.sub("\s+", '', s)
        # data_hash[name]["Seq"] = s1
        data_hash[name]["Length"] = str(len(s1))

DB_FILE.write(Redis_trans(data_hash))

os.system( "cat %s | redis-cli -n %s --pipe"%(  DB_FILE.name,db_number  ))
