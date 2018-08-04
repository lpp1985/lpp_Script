#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/4/8
'''Usage python combine_tpm   [SEQ]  [Express_appendix]''' 
from lpp import *
all_stats_file = glob.glob( '*.%s'%( sys.argv[2] ) )
all_hash = {}
all_contig = {}
total = {}
for each_f in all_stats_file:
	all_hash[ each_f.split('.')[0] ] = ''
	exec('''%s = File_dict(  open(  '%s' ,'rU'  )  ).read(1,2)'''%( each_f.split('.')[0] , each_f )  )
	for key in eval( each_f.split('.')[0] ):
		all_contig[ key ] = ''
	all_count = sum( [float(  x  )  for x in re.findall(  '\t(\S+)',open(  each_f ) .read() )] )
	total[  each_f.split('.')[0] ] = all_count
	
END = open('all_static-TPM','w')
ALL_FASTA = fasta_check(  open( sys.argv[1] ,'rU' ) )
all_leng = {}
for t,s in ALL_FASTA:
	t = t[1:-1].split()[0]

	all_leng[ t ] = len(re.sub( '\s+','',s ) )
END.write(  'Contig_id\t%s\n'%( '\t'.join( sorted(  all_hash  ) )   )  )
for each_c in sorted(all_contig):
	#if all_leng[ each_c ]<200:
		#continue
	END.write(  each_c  )
	for key in sorted( all_hash ):
		
		if each_c in eval(  key  ):
			END.write( '\t%s'%( 1000000*float(  eval( key )[  each_c ]  ) / total[ key  ]) )
		else:
			END.write(  '\t0.0'   )
	END.write('\n')