#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2016/4/22
"""
from lpp import *

def RUN_MAPPING( file_list  ):

    output_prefix = output_prefix+os.path.basename(file_list).rsplit(".",1)[0]
    bam_path = output_path+'/bam/'
    stringtie_path = output_path+'/stringtie/'
    stringtie_prefix = stringtie_path+output_prefix
    bam_prefix = bam_path+output_prefix
    read1,read2 = sorted( file_list )
    os.system(  """
     STAR --runThreadN 41 --readFilesIn %s %s --outSAMtype BAM SortedByCoordinate  --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 20 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 --twopassMode Basic  --outFileNamePrefix %s --genomeDir %s    --outSAMstrandField intronMotif  RemoveNoncanonical   && stringtie  %s*.bam  -o %s   
    
    
    """%( read1,read2,bam_prefix,refFolder,bam_prefix,stringtie_path  )  )



if __name__ == '__main__':
    '# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
    usage='''usage: python %prog [options]

    multiproecssing blast '''
    parser = OptionParser(usage =usage )

    parser.add_option("-i", "--Input", action="store",
                      dest="path",
                      type='string',
                      help="Input File")	
    
    parser.add_option("-R", "--Ref", action="store",
                      dest="refseq",
                      type='string',
                      help="refseq") 
    
    parser.add_option("-R", "--Ref", action="store",
                      dest="refseq",
                      type='string',
                      help="refseq")    
    




    parser.add_option("-F", "--Folder", action="store", 
                      dest="Folder",
                      help="RefFolder")
    
    parser.add_option("-i", "--Input", action="store", 
                      dest="InputPath",
                      help="InputPath")
    
  
    parser.add_option("-o", "--Out", action="store", 
                      dest="Out",
                      help="Out") 
    output_path = options.Out
    os.makedirs(output_path)
    
    (options, args) = parser.parse_args()
    output_path = options.Out
    refFolder = os.path.abspath( options.Folder )
    refseq = os.path.abspath( options.refseq )
    if  not os.path.exists(  options.Folder  ):
        os.makedirs(optopns.Folder)
        os.system( "STAR  --runMode genomeGenerate   --runThreadN  40 --genomeFastaFiles  %s   --genomeDir %s "%( refseq, refFolder   ) )
    output_hash = Ddict()
    for a,b,c in os.walk(  optiong.inputpath ):
        for each_f in c:
            if re.search('.pair\d+',each_f   ):
                name = each_f.split('.')[0]
                output_hash[  name ] [  a+'/'+each_f ] = ''
    out_list = []
    for each_data in output_hash:
        out_list.append( list(  output_hash[  each_data ]  )  )  
    
        
