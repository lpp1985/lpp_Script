#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/4/11
"""

from lpp import *
from optparse import  OptionParser
class Seq_All(object):
    def __init__(self,name,nul_seq,pep_seq):
        self.name = name
        self.source = self.name.split('_')[0]
        self.nul_seq = nul_seq
        self.pep_seq = pep_seq
        self.i = 0
    def __iter__( self ):
        return self
    def next ( self ):
        if self.i < len(self.pep_seq):
            location = self.i
            self.i+=1
            return self.source,self.name,location+1,self.pep_seq[location],self.nul_seq[3*location: 3*(location+1)  ]
        else:
            raise StopIteration
            
class Check_Seq(object):
    def __init__(self, nal_seq,pro_seq,END):
        self.END = END
        self.cluster =  os.path.basename( nal_seq  ).split( "." )[0] 
        
        self.all_data = []
        NAL = fasta_check( open(nal_seq,'rU')  )
        PEP = fasta_check( open(pro_seq,'rU')  )
        all_nul= {}
        all_pep = {}
        for t,s in NAL:
            s = re.sub("\s+", '', s)
            all_nul[ t.split()[0][4:]  ] = s
        for t,s in PEP:
            s = re.sub("\s+", '', s)
            all_pep[ t.split()[0][4:]  ] = s            
        for key in all_nul:
            
            self.all_data.append( Seq_All( key, all_nul[key], all_pep[key] ) )

    def __iter__(self):
        return self
    def next(self):
        return [ x.next() for x in self.all_data ]
            

    def check_data( self, all_data,all_function,all_sample):
        if self.cluster not in all_data:
            return ""
        for each_data in self:
            source_gene = Ddict()
         
            gene_aa = {}
            gene_nul = {}

            for each_loc_detail in each_data:
                source_gene[  each_loc_detail[0]  ][ each_loc_detail[1]  ]=""

                gene_aa [  each_loc_detail[1]  ]= each_loc_detail[3]
                gene_nul [  each_loc_detail[1]  ]= each_loc_detail[4]
                location = each_loc_detail[2] 
               
            if self.cluster in all_data and location in all_data[ self.cluster]:

                for sublocation in all_data[ self.cluster ][ location  ]:
                    result_data = all_data[ self.cluster ][ location  ][ sublocation  ]
                    
                    for sample in sorted( all_sample   ):
                        result_data+='\t'
                        cache_data = []
                        if sample in source_gene:
                            for gene in source_gene[sample]:
                                function = all_function[ gene ]
                                if sublocation >0:
                                    cache_data.append( gene+':'+ gene_nul [  gene  ][sublocation-1]   +' --> ' + gene_aa [  gene  ] )
                                else:
                                    cache_data.append( gene+':'+ gene_nul [  gene  ]   +' --> ' + gene_aa [  gene  ] )
                                
                        result_data += '; '.join( cache_data )
                        
                    self.END.write( result_data+'\t'+function+'\n' )
                            
                
                    
        
    

if __name__ == '__main__':


    usage = '''usage: python2.7 %prog'''
    parser = OptionParser(usage =usage ) 
    parser.add_option("-n", "--NAL", action="store", 
                      dest="NAL", 
                      default = "",
                      help="Nal Path")
    
    parser.add_option("-i", "--INPUT", action="store", 
                      dest="input_path", 
                      help="input_path")		
    parser.add_option("-o", "--Output_path", action="store", 
                      dest="output_path", 
                      help="output_path")    
    (options, args) = parser.parse_args() 
    
    
    input_path = options.input_path
    all_function ={}
    all_sample = {}
    for e_f in glob.glob(options.input_path+"/*.function"):
        FUN = open(e_f)
        sample_name = os.path.basename(e_f).split(".")[0]
        all_sample[ sample_name  ] = ""
        for line in FUN:
            line_l = line[:-1].split("\t")
            all_function[  line_l[0] ] = line_l[-1]
            
            
    OROTHO = open(options.output_path+"1.Orthologs_Cluster.txt",'rU')
    OROTHORESULT = open(options.output_path+"1.Orthologs_Cluster.txt1",'w')
    OROTHORESULT.write(OROTHO.next()[:-1]+'\tFunction\n')
    for line in OROTHO:
        line_l = line[:-1].split("\t")
        
        for data in line_l[-2:]:
            for gene in data.split(","):
                if gene in all_function:
                    function =all_function[gene]
        OROTHORESULT.write(  line[:-1]+'\t'+ function+'\n' )
    
    all_nal = glob.glob(  options.NAL+'/'+'*.nal'  )
    output_path = options.output_path
    cds_has_result = options.output_path+'/3.CDS.variation.txt'
    
    
    
    CDS_END = open( options.output_path+'/3.CDS.variation.txt2','w')
    
    RAW = open(cds_has_result,'rU')
    CDS_END.write(   RAW.next()[:-1]+'\t'+'\t'.join( sorted(   all_sample ) )+'\tFunction\n'  )
    
    all_data = Ddict()
    for line in RAW:
        line_l = line.strip().split("\t")
        cluser = line_l[0]
        location = line_l[3]
        sublocation = 0
        if "." in location:
            location,sublocation = location.split(".")
        location = int(location)
        sublocation = int( sublocation)
        all_data[ cluser ][location][sublocation] = line.strip()    
    

    
    for each_nal in all_nal:

        each_pal = each_nal.replace(".nal",".pal")
        need_check_data = Check_Seq( each_nal,each_pal,CDS_END)

        need_check_data.check_data( all_data,all_function ,all_sample)
        
