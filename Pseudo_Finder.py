#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2015/10/19
"""
from lpp import *
import os
from optparse import OptionParser
def check_path( path ):
    if not os.path.exists(path):

        os.makedirs( path )
    return os.path.abspath(path)+'/'
        
def GBLASTA( protein,assemblyresult,output ):
    #os.system("""makeblastdb  -in %s -title Assem  -parse_seqids  -out Assem -dbtype nucl"""%(assemblyresult))
    COMMAND = open("gblasta_run.bat",'w')
    RAW = fasta_check(open(protein,'rU'))
    i=0
    for t,s in RAW:
        i+=1
        
    COMMAND.write("""
    genblast -P blast  -q $input  -t %s -o $output
    """%(assemblyresult))
    os.system("""
     Genblast_Run.py -i %s -s %s  -c  %s -o %s
    """%(
       protein,COMMAND.name, i,output
       )
       ) 
def ParseGblasta(gbaresult,genewiseruncommand):
    COMMAND = open(genewiseruncommand,'w')
    cache_path = check_path("CACHE/")
    i=0
    
    data_cache_hash = {}
    GBA = block_reading(open(gbaresult,'rU'), re.escape("//******************END*******************//") )
    i=0
    for e_b in GBA:
        i+=1
        k=0
        gb_block = re.split("\n\n+", e_b)
        
        if "for query:" not in e_b:
            continue
        proteinid = re.search("for query\:\s+(\S+)", e_b).group(1)

        for align in gb_block[1:]:
            if "gene cover" not in align:
                continue
            
            
            aligndata = re.search("cover\:\d+\((\S+)\%\)\|score:([^\|]+)", align)
            perc = float(aligndata.group(1))
            score = float(aligndata.group(2))
            if perc >=80:
                i+=1
                if i not in data_cache_hash:
                    PRO= open(cache_path+'%s.pep'%(i),'w')
                    PRO.write(proteinseqHash[proteinid])
                       
                    data_cache_hash[i] = [PRO.name]
                k+=1
                NUC = open(cache_path+'%s_%s.nuc'%(i,k),'w')
                align_detail = align.split("\n")[0]
                align_detail_list = align_detail.split("|")
                subject_detail = align_detail_list[1]
                scaffold_name = subject_detail.split(":")[0]
                direct = align_detail_list[2]
                scaffoldStart,scaffoldEND = subject_detail.split(":")[1].split("..")
                scaffoldStart=int(scaffoldStart)
                scaffoldEND = int(scaffoldEND)
                if scaffoldStart<10000:
                    scaffoldStart = 0
                else:
                    scaffoldStart =scaffoldStart -10000
                scaffoldEND = scaffoldEND+10000
                
                NUC.write(">"+scaffold_name+"__%s\n"%(scaffoldStart)+assemblyseqHash[scaffold_name][scaffoldStart:scaffoldEND]+'\n')
                commandline = """Genewise_Psuedeo.py -p %s  -n %s   -o  %s.result.gff"""%(PRO.name,NUC.name,i)
                if direct =="-":
                    commandline += " -d"
                COMMAND.write(commandline+'\n')
    COMMAND.close()    
    os.system(  "cat %s | parallel -j 64"%(COMMAND.name) )  
    os.system( "cat *.result.gff > %s"%(output) )
    os.system(" rm *.result.gff")
                    

    #os.system("cat %s| parallel -j %s >genewise.out")
    
if __name__=='__main__':
    usage = '''usage: python2.7 %prog [options] Kmer




    Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
    parser = OptionParser(usage =usage )



    parser.add_option("-c", "--CPU", action="store",
                      dest="cpu",
                      type='int',
                      default = 60,
                      help="CPU number for each thread")



    parser.add_option("-p", "--pro", action="store",
                      dest="protein",
                      help="protein sequence!!")
    parser.add_option("-a", "--assembly", action="store",
                      dest="assembly",
                      help="Assemblied Genome!!")
    parser.add_option("-o", "--out", action="store",
                      dest="output",
                      default = 'genewise.out',
                      help="The output file  you want!!")
    (options, args) = parser.parse_args()
    cpu = options.cpu
    protein = options.protein
    assembly = options.assembly
    output = options.output
    
    assemblyseqHash = {}
    for t,s in fasta_check(open(assembly,'rU')):
        t  = t.split()[0][1:]
        s = re.sub("\s+",'',s)
        assemblyseqHash[t]=s
    proteinseqHash = {}
    for t,s in fasta_check(open(protein,'rU')):
        proteinseqHash[t.split()[0][1:]] = t+s
        
    GBLASTA(protein, assembly,"geneblasta.out")
    ParseGblasta("geneblasta.out", "genewise.command")    
    os.remove("genewise.command")
    os.system("rm CACHE -rf")
    os.system("rm cache -rf")
    os.system( "rm *.xml")
        
    