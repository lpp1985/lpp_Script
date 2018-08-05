#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose:
# Created: 2011/5/9

from lpp import *
from optparse import OptionParser
import os


class fastq_quality_class(  fastq_check    ):
    @staticmethod
    def check_file(  fil  ):
        ''' check the file if eixts'''
        if not os.path.isfile( fil ) and not os.path.islink( fil ):
            print(  'File %s does not exist!!'%( fil  ) )
            raise ValueError
    @staticmethod
    def check_path( total_path , check= False ):


        if check:
            '''Check the abspath, if it does not exist , It will automatically report error'''
            [ path, fil ]=os.path.split(  total_path  )
            path = os.path.abspath(  os.path.dirname(   path  )   )
            if os.path.exists( path ) and os.path.isfile( total_path ):
                return [ path, fil ]
            else:
                print( 'Path %s does not exit!!'%(path)  )
                raise ValueError
        if not os.path.exists( total_path  ):
            os.makedirs( total_path  )

        return os.path.abspath( total_path  )+os.sep

    def __init__( self,**kword  ):
        RAW  = open( os.path.abspath( kword['input_data'] ), 'rU' )
        self.RAW = super(  fastq_quality_class,self     ).__init__( RAW   )
        self.filter_read_count =0
        self.filter_data_count = 0.0
        self.raw_reads_count=0
        self.raw_data_count = 0.0
        [ root_path, file_name ]  = self.check_path( kword[ 'input_data'  ] ,check = True    )
        self.file_name = file_name
        self.sample_name  = file_name.split('.')[0]



        self.stats_path = root_path+'/stats/'
        self.check_path(  self.stats_path  )
        self.filter_path = root_path+'/Filter/%s/'%( self.sample_name )
        self.check_path(  self.filter_path )
        self.GC = 0
        #self.AT = 0
        self.N = 0
        self.Q20 = 0


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


        self.fail_file = open( self.filter_path+ file_name+ '.failFilter' ,'w')
        self.trim = False
        if  kword['trim']:
            self.trim_path = root_path+'/Trim/%s/'%( self.sample_name )


            self.check_path(  self.trim_path   )
            self.trim_data_count = 0.0
            self.trim_read_count = 0
            self.trim = True
            self.length = kword[ 'trim_length'  ]
            self.succ_trim = open( self.trim_path + file_name+'.Trim'  ,'w')
            self.fail_trim = open( self.trim_path + file_name +'.failTrim' ,'w')
            self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tTrimReads\tTrimData(MB)\tTrimReadPerc(%)\tTrimDataPerc(%)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)'
        else:
            self.stats_title = self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)'
        self.output = open( self.filter_path+  file_name +'.Filter', 'w'    )

        ## decide the threshold is locate in (0,1), then it should represent low quality percentage
        ## else it represent the total number of low quality base insequence,the thre_status is 'number'
        ## represents this
        ## situation else its value is none

    ## This function's input is  quality line of fastq and output is True of False of this reads could be reserved
    ## quality is the quality threshold to filter, quality_line in the quaility_line of fastq, low_quality is
    ## the low_quality base number of fastq,low_threshold is the threshold variant to record threshold , if
    ## thre-status == number, it represents filter threshold is low quality base number else it is low quality
    ## percentage
    def threshold_check( cls,quality_line   ):
        quality_line = re.sub( '\s+','' ,quality_line)
        total_length = len( quality_line )
        low_quality = len(  [ x     for x in quality_line if cls.phred_64_quality_hash[x]<cls.quality    ] ) ## this is the base number of  total low quality base

        if cls.threshold > 1:
            low_threshold = low_quality
        else:
            low_perc = float( low_quality  ) / total_length
            low_threshold = low_perc
        if low_threshold > cls.threshold:
            return False
        else:

            return True
        ## This function plot the quality picture of each run ##
    @staticmethod
    def quality_plot( FILE_HANDLE    ):

        os.system(  '''( fastx_quality_stats -i %s > %s.qualstats &&'''%(  FILE_HANDLE.name,FILE_HANDLE.name  ) + '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png&&'''%(  FILE_HANDLE.name,FILE_HANDLE.name   ) + '''fastx_nucleotide_distributionPer.R  %s.qualstats %s.nucdistr.png ) &'''%(  FILE_HANDLE.name,FILE_HANDLE.name    ) )




    def mutiprocess_plot(  self  ):


        all_file_attribure = []
        for each_attrib in dir(self):

            attribure_name  = eval(  'self.'+ each_attrib  )
            if isinstance( attribure_name,file  ):
                all_file_attribure.append(  attribure_name  )


        map( self.quality_plot , all_file_attribure )



    def next(self):
        ## The Try excpt capture the last iteration of programe and output the quality of each_file
        ## output_hash is a variant record reads_status of a read,status is variant to record reads status,the value
        ## name is the reads name. The value failTrim is trim fail-status, the value Trim is represents could cold Trim
        ## the value fail_Filter is low quality Filter , value Filter is could be successfuly filter Reads

        output_hash = {}

        status = ''


        try:
            name,seq,define2,quality = super( fastq_quality_class, self  ).next(   )

            raw_cont =''.join( [  name,seq,define2,quality  ]   )

            self.GC += len(    re.findall( '([G|C])',seq    )     )

            self.N += len(  re.findall( '(N)',seq    )    )

            self.Q20 += len(  [ x     for x in re.sub( '\s+$','',  quality ) if self.phred_64_quality_hash[x]>=20    ] )

            self.raw_data_count+= len( seq[:-1]  )
            self.raw_reads_count +=1

        except StopIteration:
            self.stats = {}

            STATS = open( self.stats_path + re.sub( '(\.[^\.]+)$','', self.file_name   ) +'.stats'  ,'w')
            [name,category ]= self.file_name.split('.')[:2]



            output_rank = 	[ '''\'%s\'%( self.raw_reads_count )''',
                            '\'%.2f\'%( float(  self.raw_data_count  ) /1048576  )',
                            '\'%s\'%( self.trim_read_count )',
                            '\'%.2f\'%( float(  self.trim_data_count )/1048576  )',
                            '\'%.2f\'%(100 * self.trim_read_count/float( self.raw_reads_count )  )',
                            '\'%.2f\'%( 100 * self.trim_data_count/self.raw_data_count, )',
                            '\'%.2f\'%( self.filter_read_count )',
                            '\'%.2f\'%( float(  self.filter_data_count )/1048576 )',
                            '\'%.2f\'%( 100 * self.filter_read_count/float( self.raw_reads_count   ) )' ,
                            '\'%.2f\'%( 100*self.filter_data_count/float(   self.raw_data_count   ) ) ' ,
                            '\'%.2f\'%( 100*self.N/float(   self.raw_data_count   ) ) ' ,
                            '\'%.2f\'%( 100*self.GC/float(   self.AT + self.GC   ) ) ' ,
                            '\'%.2f\'%( 100*self.Q20/float(   self.raw_data_count   ) ) '

                            ]

            output = []
            for each_data in output_rank:
                try:

                    output.append( eval( each_data )  )
                except:
                    pass
            output_data = '\t'.join( output   )

            STATS.write(  self.stats_title+'\n'  )
            STATS.write( name+'\t'+category+'\t'+output_data+'\n' )

            self.mutiprocess_plot()

            raise  StopIteration

        ### Trim the reads by the N Reads locate in  Head or  Tail of
        output_hash[ 'name' ] = name[:-3]
        if self.trim:
            seq_re = re.search( '^N*(\S+?)N*\n$' ,seq     )
            if  seq_re:
                seq_locat     = seq_re.span(1)
                seq_trim      = seq[ seq_locat[0]:seq_locat[1]  ]
                quality_trim  = quality[ seq_locat[0]:seq_locat[1  ] ]
                output_cont   = name+seq_trim+'\n'+define2+quality_trim+'\n'



                if len(  seq_trim   ) >= self.length:
                    raw_cont            = output_cont
                    seq                 = seq_trim
                    quality             = quality_trim

                    self.trim_data_count += len( seq_trim )
                    self.trim_read_count += 1

                    self.succ_trim.write(  output_cont   )

                    status = 'Trim'

                    output_hash[ status ] = output_cont




                else:

                    self.fail_trim.write(  raw_cont   )

                    status = 'failTrim'

                    output_hash[ status ] = raw_cont


                    return (  output_hash  )

            else:
                status = 'failTrim'

                output_hash[ status ] = raw_cont

                self.fail_trim.write(  raw_cont   )

                return (   output_hash     )


        if self.threshold_check(  quality   ):
            self.filter_read_count +=1
            self.filter_data_count += len( seq[:-1] )
            self.output.write(  raw_cont  )
            status = 'Filter'
            output_hash[ status ] = raw_cont
        else:
            status = 'failFilter'
            output_hash[ status ] = raw_cont
            self.fail_file.write(  raw_cont   )

        return output_hash

#if __name__ =='__main__':


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

all_raw_file = glob.glob( path+'*.'+ options.appendix)







trim = options.trim
if trim :
    status = trim
    trim_length = options.trim_length

quality = options.quality
threshold = options.threshold
if __name__ =='__main__':
    for each_f in all_raw_file:
        input_data = each_f
        fastq = fastq_quality_class( input_data = input_data,trim = trim, quality = quality, threshold = threshold ,trim_length = trim_length )
        for x in fastq:
            pass