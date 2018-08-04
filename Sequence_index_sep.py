#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/9/7
from lpp import *
from optparse import OptionParser
#from multiprocessing import Manager

# you could type [  SCRIPT_NAME  ] -h to see the help log !!!!
usage='''usage: python %prog [options] '''
parser = OptionParser(usage =usage )
parser.add_option("-p", "--Path", action="store", 
                  dest="path",
                  type='string',  
                  default = './',
                  help="Input path")
parser.add_option("-i", "--Index", action="store", 
                  dest="index",
                  type='string',  
                  help="Index Fasta File")
parser.add_option("-s", "--Score", action="store", 
                  dest="score",
                  type='int',  
                  default = 15,
                  help="Bit Score to recognize Reads")
parser.add_option("-t", "--Tolerate", action="store", 
                  dest="tolerate",
                  type='int',  
                  default  = 20,
                  help="Number to consider the index is too similar to recognize")


parser.add_option("-o", "--Output", action="store", 
                  dest="output",
                  default  = './',
                  help="Output Path")

(options, args) = parser.parse_args() 
#输入数据的目录
path = options.path
if path[-1]!='/':
    path +='/'
#    Index所在文件，要求是Fasta格式的
#    如 >Index1
#       ATCGGG
index = options.index


# score 是代表对序列进行分类的时候认为序列相似的最小阈值，小于这个值，即使是best hit one，也认为是不相似的
score = options.score

# tolerate 是代表index之间相似度的衡量标准，如果过于相似，就会导致序列无法分开
# 因此，设置一个容忍值，如果过于相似，改程序将拒绝运行，代表该数据分不开

tolerate = options.tolerate



#output 代表输出目录

output = options.output




########################################################################################

#INDEX为文件句柄，指打开的fasta格式的Index文件

INDEX = fasta_check(  open(  index,'rU'  ) )

index_hash = {}
all_index = []
for t,s in INDEX:
    s = re.sub( '\s+','',s )
    index_hash[s] = t[1:-1]
    all_index.append( s )
    
#衡量index是否过于相似，导致无法分离
for each_seq1 in all_index:
    number = all_index.index( each_seq1 )
    for each_seq2 in all_index[ number+1:  ]:
        end_score = SmithWaterman( each_seq1,each_seq2  ).align()

        if end_score >=tolerate:
            print( 
                'Error！！！%s and %s is too same to seperate'% (index_hash[  each_seq1 ] , index_hash[  each_seq2 ]  ) 
            )
            sys.exit()


#找到所有的数据并分类
#File_Hash是记录文件结果
File_Hash = {}

all_file = glob.glob(  path+'*.out'  )

for each_f in all_file:
    Tag = re.search( 'P(\d+)_',each_f   ).group(1)
    if Tag == '1':
        Categary = "Read1"
    elif Tag=='3':
        Categary = "Read2"
    else:
        Categary = "Index"
        INDEX_FILE = fastq_check(  open( each_f,'rU' )  )
    File_Hash[ Categary  ] =fastq_check(  open( each_f,'rU' )  )
    
class INDEX_SEP( object ):
    def __init__( self,File_Hash,Index_hash, Output   ):
        Index_hash[ '' ] = 'Error'
        self.Index_hash = Index_hash
        #Output为输出目录
        if Output[-1] !='/':
            Output+='/'
        #output_file_hash用于存储输出目录，为{ 样本：   { ‘Read1’：输出目录  ，‘Read2’：数据处目录 }}
        output_file_hash = Ddict()
        for each_catergory  in File_Hash:

            setattr(  self,each_catergory,File_Hash[  each_catergory ]  )
            for seq, name in index_hash.items():
                output_file_hash[  name  ][  each_catergory ] = open( Output+name+'.%s.fastq'%( each_catergory   )  ,'w')
        
        self.output_file_hash = output_file_hash

    def get_allIndex( self ):
        #将所有的index出现的序列进行统计
        all_index_hash = Ddict()
        for a,b,c,d in INDEX_FILE:
            seq = b
            seq_name = a[:-2]
            all_index_hash[ seq ] [ seq_name ] = ''
        self.all_index_hash = all_index_hash
       

        
    def classifer(  self,query_seq  ):
        '''
        对index进行分类的程序,用于对单个index进行分类，如果序列打分超过一定的阈值，则认为该序列没问题，然后根据排名查看两个index打分是否相同，如果相同或者打分都很低，则认为该Index有问题，分不出来！否则，输出分类结果
        '''
        score_hash = Ddict()
        for each_index_seq in self.Index_hash:
            if each_index_seq =='':
                continue
            align_score = SmithWaterman(  query_seq,  each_index_seq  ).align()
            score_hash[ align_score  ][  self.Index_hash[ each_index_seq  ]  ] = ''
        max_score = max( score_hash )

        if max_score < score:
            return 'Error'
        elif len( score_hash[  max_score   ] )>1:
            return 'Error'
        else:
            return score_hash[ max_score  ].keys()[0]
    def run( self  ):
        '''开始分类'''
        self.get_allIndex()
        self.belong = {}
        for each_index in self.all_index_hash:
            
            belong_to = self.classifer(  each_index )
            for each_reads in self.all_index_hash[   each_index ]:
                
                self.belong[  each_reads   ]= belong_to
        for a,b,c,d in self.Read1:
            name = a[:-2]
            write_category = self.belong[ name  ]
            self.output_file_hash[ write_category  ][ 'Read1' ].write( a+b+c+d  )
            
            a2,b2,c2,d2 = self.Read2.next()
            self.output_file_hash[ write_category  ][  'Read2' ].write( a2+b2+c2+d2  )
            
            ai,bi,ci,di = self.Index.next()
            self.output_file_hash[ write_category  ][  'Index' ].write( ai+bi+ci+di  )

need_sep = INDEX_SEP( File_Hash,index_hash, output )
need_sep.run()