#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/11/22

from lpp import *
from optparse import OptionParser

if __name__=='__main__':
    usage = '''usage: python2.7 %prog [options] Kmer




    Kmer is a list of K value you want,e.g  [ 1, 2, 3, 4 ]'''
    parser = OptionParser(usage =usage )



    parser.add_option("-m", "--MAF", action="store",
                      dest="maf",
                      help="maf file")

    parser.add_option("-v", "--VCF", action="store",
                      dest="vcf",
                      help="Vcf file")
    (options, args) = parser.parse_args()
    
    
    MAF =  block_reading( open( options.maf,'rU'  ),tag = 'a score='  )
    VCF = open(  options.vcf ,'w'     )
   
    VCF.write( '''##fileformat=VCFv4.1
##fileDate=20130702
##reference=NC_000913
##INFO=DP,1,Integer,"Total Depth of Coverage"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
''' )




    MAF.next()


    for block in MAF:
        line_list = block.split('\n')
        target_line, query_line = line_list[0], line_list[1]
        target_list = target_line.strip().split()
        query_list = query_line.strip().split()
        query_name,query_start,query_direction,query_alignment = query_list[1],query_list[2],query_list[4],query_list[-1]
        target_name,target_start,target_direction,target_alignment = target_list[1],target_list[2],target_list[4],target_list[-1]
        query_alignment = iter( query_alignment.upper() )
        query_start = int( query_start )
        target_start = int( target_start )
        query_tag = ''
        subject_tag  = ""
        locus = target_start
        status = 0
        target_alignment = target_alignment.upper()
        query_alignment = query_alignment
        for key in target_alignment:
            sub_key = query_alignment.next()
            if key == sub_key:
                if status:
                    VCF.write( 
                        '\t'.join(  [target_name,'%s'%( locus  ),'.',query_tag,subject_tag,'100.00\tPASS\tDP=100'   ] )+'\n'
                        
                    )
                    status = 0
                
                target_start+=1
                locus = target_start
                query_tag = key
                subject_tag = key
                
                
            else:
                status = 1

                if re.match( '\w',key ): 
                    target_start+=1
                    query_tag +=key
                if re.match( '\w',sub_key  ):
                    subject_tag +=sub_key
        else:

            if status:
                VCF.write( 
                    '\t'.join(  [target_name,'%s'%( locus  ),'.',query_tag,subject_tag,'100.00\tPASS\tDP=100'   ] )+'\n'
                    
                )
                status = 0            
