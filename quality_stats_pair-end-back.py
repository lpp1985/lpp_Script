#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/8/29
from lpp import *
import multiprocessing
from optparse import OptionParser
'# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!'
usage='''usage: python %prog [options] 

It can stats quality and filter bad quality reads'''
parser = OptionParser(usage =usage ) 
parser.add_option("-p", "--Path", action="store", 
                  dest="path",
                  type='string',  
                  help="Input path")
parser.add_option("-a", "--Appendix", action="store", 
                  dest="appendix",
                  type='string',  
                  help="Input FASTQ")

parser.add_option("-o", "--OUTPUT", action="store", 
                  dest="output", 
                  help="The fastq file successfully filter the quality")





parser.add_option("-t", "--Trim", action="store_true", default= False,
                  dest="trim", 
                  help="Do you want to trim the fastq file( trim the N in head or in tail  )")


parser.add_option("-l", "--Trim_length", action="store", 
                  dest="trim_length", 
                  type= "int",
                  help="The length of threshold  to  filter after trim, THIS PARAMATER OUGHT TO BE INGORANCE IN [-t] IS MISSING!!!")


parser.add_option("-q", "--Quality", action="store", type = "int",
                  dest="quality", 
                  help="Minimum quality score to keep.")


parser.add_option("-r", "--Threshold", action="store", type = "float",
                  dest="threshold", 
                  help='''The Threashold you want to filter, if this paramater's value locate in the area of (0,1), It represent the maxium percentage the low quality accounts for , 
                  if Integer it represents the maxium number of low quality base in the whole reads to trim''')

(options, args) = parser.parse_args() 

path = options.path

## Split file into different group

all_raw_file = glob.glob( path+'/*.'+ options.appendix)

all_file_group = {}

for each_f in all_raw_file:
    name = re.search( '^([^\.]+)',os.path.split(os.path.abspath(each_f)  )[-1] ).group(1)
    name_suffix = re.escape( name )
    all_file_group[ name] = [  x   for x in all_raw_file if re.search( '^('+name_suffix+')',  os.path.split(os.path.abspath(x)  )[-1]  )     ]





trim = options.trim
if trim :
    status = trim
    trim_length = options.trim_length

quality = options.quality
threshold = options.threshold



class fastq_quality_class(  fastq_check    ):
    @staticmethod
    def check_path( path ):
        if not  os.path.exists(path):

            os.makedirs( path )
    def __init__( self,kword  ):
        abs_path = os.path.abspath( kword['path'] + kword[ 'input_data'  ]    )
        [root_path, file_name ]  = os.path.split( abs_path  )
        sample_name  = file_name.split('.')[0]

        trim_path = root_path+'/Trim/%s/'%( sample_name )

        filter_path = root_path+'/Filter/%s/'%( sample_name )

        self.check_path(  trim_path )
        self.check_path(  filter_path )

        RAW  = open( kword['input_data'] , 'rU' )
        super(  fastq_quality_class,self     ).__init__( RAW   )

        ## the quality hash to represent the quality in each site,which can make searching more fast
        self.phred_64_quality_hash={}
        for x in xrange(1,41):
            self.phred_64_quality_hash[ chr( x+64  ) ] = x
        ## self.attribute record the variant input,quality ->quality threshold, threshold -> to threshold you 
        ##want, length -> reads_length after trim to recieve, fail_file -> low quality reads, trim-> True: to 
        ##trim, False: Not trim, succ_trim -> trim successfuly file to record, fail_trim -> trim failed reads to 
        ##record, output-> file which filtered reads write in 
        self.quality = kword['quality']
        self.threshold = kword['threshold']

        self.fail_file = open( filter_path+kword[ 'input_data'  ] +'.failFilter' ,'w')
        self.trim = False
        if  kword['trim']:
            self.trim = True
            self.length = kword[ 'trim_length'  ]
            self.succ_trim = open( trim_path+kword[ 'input_data'  ] +'.Trim'  ,'w')
            self.fail_trim = open( trim_path+kword[ 'input_data'  ] +'.failTrim' ,'w')
        self.output = open( filter_path+  kword[ 'input_data'  ] +'.Filter', 'w'    )

        ## decide the threshold is locate in (0,1), then it should represent low quality percentage
        ## else it represent the total number of low quality base insequence,the thre_status is 'number' 
        ## represents this 
        ## situation else its value is none
        if threshold >1:
            self.thre_status = 'number'
        else:
            self.thre_status = ''
    ## This function's input is  quality line of fastq and output is True of False of this reads could be reserved 
    ## quality is the quality threshold to filter, quality_line in the quaility_line of fastq, low_quality is 
    ## the low_quality base number of fastq,low_threshold is the threshold variant to record threshold , if 
    ## thre-status == number, it represents filter threshold is low quality base number else it is low quality 
    ## percentage
    def threshold_check( cls,quality_line   ):
        quality_line = re.sub( '\s+','' ,quality_line)
        total_length = len( quality_line )
        low_quality = len(  [ x     for x in quality_line if cls.phred_64_quality_hash[x]<cls.quality    ] ) ## this is the base number of  total low quality base
        if cls.thre_status:
            low_threshold = low_quality
        else:
            low_perc = float( low_quality  ) / total_length
            low_threshold = low_perc
        if low_threshold>cls.threshold:
            return False
        else:

            return True
        ## This function plot the quality picture of each run ##
    @staticmethod
    def quality_plot( FILE_HANDLE    ):

        os.system(  '''fastx_quality_stats -i %s > %s.qualstats'''%(  FILE_HANDLE.name,FILE_HANDLE.name  )   )
        os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png'''%(  FILE_HANDLE.name,FILE_HANDLE.name   )   )
        os.system('''fastx_nucleotide_distribution.R  %s.qualstats %s.nucdistr.png'''%(  FILE_HANDLE.name,FILE_HANDLE.name    ))
    

         
    def mutiprocess_plot(  self  ):
        pool = multiprocessing.Pool( processes=20 )

        all_file_attribure = []
        for each_attrib in dir(self):

            attribure_name  = eval(  'self.'+ each_attrib  )
            if isinstance( attribure_name,file  ):
                all_file_attribure.append(  attribure_name  )

        
        map( self.quality_plot , all_file_attribure )
        
        
        
    def next(self):
        ## The Try excpt capture the last iteration of programe and output the quality of each_file
        ## output_hash is a variant record reads_status of a read,status is variant to record reads status,the value
        ##name is the reads name. The value failTrim is trim fail-status, the value Trim is represents could cold Trim
        ## the value fail_Filter is low quality Filter , value Filter is could be successfuly filter Reads 

        output_hash = {}

        status = ''

        try:

            name,seq,define2,quality = super( fastq_quality_class, self  ).next(   )

        except:

            self.mutiprocess_plot()

            raise  StopIteration

        ## Trim the reads by the N Reads locate in  Head or  Tail of
        output_hash[ 'name' ] = name[:-3]
        if self.trim:
            seq_re = re.search( '^N*(\S+?)N*\n$' ,seq     )
            if  seq_re:
                seq_locat = seq_re.span(1)
                seq_trim = seq[ seq_locat[0]:seq_locat[1]  ]
                if len(  seq_trim   ) >= self.length:

                    seq = seq_trim+'\n'

                    quality = quality[ seq_locat[0]:seq_locat[1  ] ] +'\n'

                    self.succ_trim.write(  name+seq+define2+quality   )

                    status = 'Trim'

                    output_hash[ status ] = name + seq + define2 + quality
                    if self.threshold_check(  quality   ):

                        self.output.write(  name+seq+define2+quality   )
                        status = 'Filter'
                        output_hash[ status ] = name + seq + define2 + quality 
                    else:
                        status = 'failFilter'
                        output_hash[ status ] = name + seq + define2 + quality 
                        self.fail_file.write(  name+seq+define2+quality   )


                else:

                    self.fail_trim.write(  name+seq+define2+quality    )

                    status = 'failTrim'

                    output_hash[ status ] = name + seq + define2 + quality 

                    
                    return (  output_hash  )

            else:
                status = 'failTrim'
                output_hash[ status ] = name + seq + define2 + quality 
                self.fail_trim.write(  name+seq+define2+quality    )
               
                return (   output_hash     )


        
        return output_hash
class fastq_quality_pair(  object    ):

    ## input a hash of a copy of  locals----input_hash = locals().copy(), the set the files variant to be your input file_list, just like input_hash[ 'files'  ] = [ read1,read2     ]
    def __init__(  self,kword  ):
        all_file = sorted(kword[ 'files' ], key= lambda x : int( re.search( 'read(\d)',x,re.I ).group(1) )      )
        read1_kword = kword.copy()
        read2_kword = kword.copy()
        read1_kword[ 'input_data'  ] = all_file[0]

        read2_kword[ 'input_data'  ] = all_file[-1]

        self.read1 = fastq_quality_class(  read1_kword  )
        self.read2 = fastq_quality_class(  read2_kword  )
        # this function change the name of raw data name to be paired-name
        def trim_name( char ):
            char = re.sub( '\.read\d.+?$','',char  )
            return char

        #add new name of paired data file
        if trim:
            self.read1_trim_pair_output   = open(   trim_name( self.read1.succ_trim.name )+'.Trim.pair1' ,'w' )
            self.read2_trim_pair_output   = open(   trim_name( self.read2.succ_trim.name )+'.Trim.pair2' ,'w'  )
        self.read1_filter_pair_output = open(   trim_name( self.read1.output.name )+'.Filter.pair1' ,'w'  )
        self.read2_filter_pair_output = open(   trim_name( self.read2.output.name )+'.Filter.pair2' ,'w'  )
        self.log = open(  trim_name ( self.read1.file_handle.name )+'.list','w'    )


    def __iter__(self):

        return self

    def next( self  ):
        def check( output_hash ):
            for key in [ 'Filter','failFilter','Trim','failTrim'    ]:
                if key in output_hash:
                    status1 = key
                    return status1
        ## the status_transfer_hash can make output more accuracy , and output_hash record the output file to be output
        status_transfer_hash = { 'Filter':'Trim' ,'failFilter' :'Trim','Trim':'Trim','failTrim':'failTrim'  }
        if trim:
            File_hash = { 'Filter':[self.read1_filter_pair_output,self.read2_filter_pair_output],'Trim':[ self.read1_trim_pair_output, self.read2_trim_pair_output      ] } 
        else:
            File_hash = { 'Filter':[self.read1_filter_pair_output,self.read2_filter_pair_output], } 
        try:
            output_hash1 = self.read1.next()
            output_hash2 = self.read2.next()
        except:
            output_hash2 = self.read2.next()

        self.log.write(  output_hash1[ 'name'   ]  )

        status1 = check(  output_hash1  )
        
        status2 = check(  output_hash2  )

        if status1== status2:
            if status1 in File_hash:
                File_hash[ status1  ][0].write( output_hash1[  status1 ]   )
                File_hash[ status1  ][-1].write( output_hash2[  status2 ]   )
        #else:
            status_transfer1 = status_transfer_hash[ status1 ]

            status_transfer2 = status_transfer_hash[ status2 ]

            if status_transfer1==status_transfer2:

                if status_transfer1 in File_hash:
                    File_hash[ status_transfer1  ][0].write( output_hash1[  status_transfer1 ]   )
                    File_hash[ status_transfer1  ][-1].write( output_hash2[  status_transfer2 ]   )
        self.log.write( '\t'+status1+'\t'+status2+'\n'   )





for each_f in all_file_group:
    files = all_file_group[ each_f ]
    aa = fastq_quality_pair( locals() )
    for i in aa:
        pass