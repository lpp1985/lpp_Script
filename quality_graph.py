#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/9/7
from lpp import *
from optparse import OptionParser
#from quality_stats_singleEnd import fastq_quality_class
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
usage='''usage: python %prog [options] '''
parser = OptionParser(usage =usage )
parser.add_option("-p", "--Path", action="store", 
                  dest="path",
                  type='string',  
                  help="Input path")
#parser.add_option("-a", "--Appendix", action="store", 
                  #dest="appendix",
                  #type='string',  
                  #help="Input FASTQ")

(options, args) = parser.parse_args() 
path = options.path
#appendix = re.escape(  options.appendix )
def quality_plot( FILE_HANDLE    ):

    os.system(  '''(fastx_quality_stats -i %s > %s.qualstats&&'''%(  FILE_HANDLE.name,FILE_HANDLE.name  ) + '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png&&'''%(  FILE_HANDLE.name,FILE_HANDLE.name   ) + '''fastx_nucleotide_distributionPer.R  %s.qualstats %s.nucdistr.png)&'''%(  FILE_HANDLE.name,FILE_HANDLE.name    ) )
    #os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  FILE_HANDLE.name,FILE_HANDLE.name   )   )
    #os.system('''fastx_nucleotide_distributionPer.R  %s.qualstats %s.nucdistr.png'''%(  FILE_HANDLE.name,FILE_HANDLE.name    ))





for a,b,c in os.walk( path ):
    for each_f in c:

        RAW = open(  a+'/'+each_f ,'rU' )

        for line in RAW:
            if line .startswith('@'):



                quality_plot(  RAW )
            break
            