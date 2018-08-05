#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/5/25
from lpp import *
FASTA = fasta_check( open( sys.argv[2],'rU'  )  )
all_seq = {}
for t,s in FASTA:
	title = t[1:-1]
	s= re.sub( '\s+','' ,s)
	all_seq[ title ] = s
END = open( 'Chromsome.fasta','w' )
end_seq = []
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
			end_seq.append( seq )
end_seq = ''.join(  end_seq )
end_seq = re.sub( '(\w{60})','\\1\n',end_seq )
END.write(  '>2M_seq\n%s\n'%( end_seq )  )
		