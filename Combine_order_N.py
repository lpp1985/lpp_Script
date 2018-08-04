#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/25
from lpp import *
FASTA = fasta_check( open( sys.argv[2],'rU'  )  )
all_seq = {}
len_1 = []
for t,s in FASTA:
	title = t[1:-1]
	s= re.sub( '\s+','' ,s)
	all_seq[ title ] = s
	len_1.append( len(s) )
print( sum( len_1 ) )
END = open( 'Chromsome_N.fasta','w' )
end_seq = []
all_order = []

len_2 = []
for block in re.split( '\n{2,}' , open( sys.argv[1],'rU' ).read()  ):
	
	if block.startswith( 'Ordered Contigs' ):
		line_all = block.split('\n')
		line_all = line_all[2:] 
		
		for line in line_all:
			
			line_l = line.split()
			name = line_l[1]
			tag = line_l[3]
			seq = all_seq[ name ]
			if tag == 'complement':

				seq = complement( seq )
			len_2.append( len( seq  ) )
			
			end_seq.append( seq )
			all_order.append( name+' '+ tag  )
print( sum(len_2) )
ORDER = open( 'Scaffold_N.order','w'  )
i=0
for each_tag  in all_order[1:]:
	i+=1
	ORDER.write( 'GAP_%s'%(i)+'\t'+all_order[ all_order.index(each_tag)-1]+'\t'+each_tag+'\n'  )
end_seq = '%s'%( 'N'*1000 ).join(  end_seq )
end_seq = re.sub( '(\w{60})','\\1\n',end_seq )
END.write(  '>2M_seq\n%s\n'%( end_seq )  )
		