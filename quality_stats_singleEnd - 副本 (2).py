#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/9
from lpp import *

from lpp import *
import multiprocessing
from optparse import OptionParser



class fastq_quality_class(  fastq_check    ):
    @staticmethod
    def check_path( path ):
        if not  os.path.exists(path):

            os.makedirs( path )
            
    def __init__( self,kword  ):
        RAW  = open( kword['input_data'] , 'rU' )
        super(  fastq_quality_class,self     ).__init__( RAW   )
        self.filter_read_count =0
        self.filter_data_count = 0.0
        self.raw_reads_count=0
        self.raw_data_count = 0.0
        abs_path = os.path.abspath( kword['path'] + kword[ 'input_data'  ]    )
        [root_path, file_name ]  = os.path.split( abs_path  )
        sample_name  = file_name.split('.')[0]

        trim_path = root_path+'/Trim/%s/'%( sample_name )

        filter_path = root_path+'/Filter/%s/'%( sample_name )
        
        self.stats_path = os.path.split(  os.path.abspath(  self.file_handle.name   ) )[0]+'/stats/'
        self.check_path(  self.stats_path  )
        self.check_path(  trim_path   )
        self.check_path(  filter_path )
        self.GC = 0
        self.AT = 0
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
        
        
        self.fail_file = open( filter_path+kword[ 'input_data'  ] +'.failFilter' ,'w')
        self.trim = False
        if  kword['trim']:
            self.trim_data_count = 0.0
            self.trim_read_count = 0
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

        os.system(  '''( fastx_quality_stats -i %s > %s.qualstats &&'''%(  FILE_HANDLE.name,FILE_HANDLE.name  ) + '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png&&'''%(  FILE_HANDLE.name,FILE_HANDLE.name   ) + '''fastx_nucleotide_distributionPer.R  %s.qualstats %s.nucdistr.png ) &'''%(  FILE_HANDLE.name,FILE_HANDLE.name    ) )
        #os.system(  '''fastx_quality_boxplot.R %s.qualstats %s.qualstats.png&'''%(  FILE_HANDLE.name,FILE_HANDLE.name   )   )
        #os.system('''fastx_nucleotide_distributionPer.R  %s.qualstats %s.nucdistr.png&'''%(  FILE_HANDLE.name,FILE_HANDLE.name    ))
    

         
    def mutiprocess_plot(  self  ):
        #pool2 = multiprocessing.Pool( processes=20 )

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
    
            self.AT +=  len(   re.findall( '([A|T])',seq    )   )
            
            self.GC += len(    re.findall( '([G|C])',seq    )     )
    
            self.N += len(  re.findall( '(N)',seq    )    )
            
            self.Q20 += len(  [ x     for x in re.sub( '\s+$','',  quality ) if self.phred_64_quality_hash[x]>=20    ] )
          
            self.raw_data_count+= len( seq[:-1]  )
            self.raw_reads_count +=1

        except:
            self.stats = {}
            self.stats_title = ''
            
            STATS = open( self.stats_path+re.sub( '(\.[^\.]+)$','', self.file_handle.name    ) +'.stats'  ,'w')
            [name,category ]= os.path.split(self.file_handle.name)[-1].split('.')[:2]
            
            
                        
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
                
   
            if self.trim :
                self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tTrimReads\tTrimData(MB)\tTrimReadPerc(%)\tTrimDataPerc(%)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)' 
            else:
                self.stats_title = '#Sample\tCategory\tRawReads\tRawData(MB)\tFilterReads\tFilterData(MB)\tFilterReadsPerc(%)\tFilterDataPerc(%)\tN(%)\tGC(%)\tQ20(%)' 
            STATS.write(  self.stats_title+'\n'  )
            STATS.write( name+'\t'+category+'\t'+output_data+'\n' )

            self.mutiprocess_plot()

            raise  StopIteration

        ### Trim the reads by the N Reads locate in  Head or  Tail of
        output_hash[ 'name' ] = name[:-3]
        if self.trim:
            seq_re = re.search( '^N*(\S+?)N*\n$' ,seq     )
            if  seq_re:
                seq_locat = seq_re.span(1)
                seq_trim = seq[ seq_locat[0]:seq_locat[1]  ]
                quality_trim = quality[ seq_locat[0]:seq_locat[1  ] ] 
                quality_end = quality_trim
                seq_end = seq_trim
                # Trim quality from  start
                for each_quality in quality_trim :
                    if self.phred_64_quality_hash[each_quality]< self.quality:
                        seq_end= seq_end[1:]
                        quality_end = quality_end[1:]
                    else:
                        break
                # Trim quality from end
                if seq_end:
                    for each_quality in quality_trim[::-1]:
                        if self.phred_64_quality_hash[each_quality]< self.quality:
                            seq_end= seq_end[:-1]
                            quality_end = quality_end[:-1]
                        else:
                            break
                    
                if len(  seq_end   ) >= self.length:
                    self.trim_data_count += len( seq_end )
                    self.trim_read_count += 1


                    self.succ_trim.write(  name+seq_end+'\n'+define2+quality_end+'\n'   )

                    status = 'Trim'

                    output_hash[ status ] = name+seq_end+'\n'+define2+quality_end+'\n'
                    if self.threshold_check(  quality_end+'\n'   ):
                        self.filter_read_count +=1
                        self.filter_data_count += len( seq_end )
                        self.output.write(  name+seq_end+'\n'+define2+quality_end+'\n'  )
                        status = 'Filter'
                        output_hash[ status ] = name+seq_end+'\n'+define2+quality_end+'\n'
                    else:
                        status = 'failFilter'
                        output_hash[ status ] = name+seq_end+'\n'+define2+quality_end+'\n'
                        self.fail_file.write(  name+seq_end+'\n'+define2+quality_end+'\n'   )


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

        else:
            if self.threshold_check(  quality   ):
                self.filter_read_count +=1
                self.filter_data_count += len( seq[:-1] )
                self.output.write(  name+seq+define2+quality   )
                status = 'Filter'
                output_hash[ status ] = name + seq + define2 + quality 
            else:
                status = 'failFilter'
                output_hash[ status ] = name + seq + define2 + quality 
                self.fail_file.write(  name+seq+define2+quality   )
        
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
        fastq = fastq_quality_class( locals() )
        for x in fastq:
            pass