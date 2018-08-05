#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: Automaticly produce report of QC
# Created: 2011/5/10
from lpp import *
from optparse import OptionParser
usage = '''usage: python %prog [options] 

It can automaticly produce QC report by WORD locate in C:/Report/'''
parser = OptionParser(usage =usage ) 

import win32com, win32com.client, sys,os,os.path
word=win32com.client.Dispatch("Word.Application")
word.Visible = True
dir_output = os.path.abspath(  'c:/'  )+'Report/'
parser.add_option("-p", "--Path", action="store", 
                  dest="path",
                  type='string',  
                  help="the path of RAW data")
(options, args) = parser.parse_args() 
inputpath = options.path
check_path( dir_output )
output_hash = Ddict()
for r,d,f in os.walk( inputpath+'/' ):
	for each_f in f:
		if '.png' in each_f:
			name  = each_f.split('.')[0]

			path = os.path.abspath( r )+'/'
			other = each_f.split('.',2)[-1]
			tag = each_f.split('.',)[-2]
			if 'other' in each_f:
				output_hash[ name ][ 'filter' ][ tag  ]= path +each_f
			else:
				output_hash[ name ][ 'raw' ][ tag  ]= path +each_f

		if '.other.stats' in each_f :

			RAW = open( path+each_f,'rU' )

			output_hash[ name ]['filterreads_content'] = RAW.next()[:-1]

			output_hash[ name ]['filterreads_number'] = RAW.next()[:-1]

		elif 'fastq.stats' in each_f:
			RAW = open( path+each_f,'rU' )
			output_hash[ name ]['raw_content'] = RAW.next()[:-1]
			output_hash[ name ]['raw_number'] = RAW.next()[:-1]



#print(  output_hash   )
print( output_hash )
def MS_Word_Find_Past_Image(app, Search_Word, fileName):
	wdShapeCenter = -999995
	wdWrapSquare = 0
	wdShapeInside                 =-999994    # from enum WdShapePosition
	wdWrapInline                  =0x7        # from enum WdWrapType
	wdCharacter                   =0x1        # from enum WdUnits

	find = app.Selection.Find
	find.Text = Search_Word
	while app.Selection.Find.Execute():
		print "I fonud ", Search_Word

		app.Selection.Delete(Unit = wdCharacter, Count = 1)

		Anchor = app.Selection.Range 
		objShape = app.ActiveDocument.InlineShapes 
		handle = objShape.AddPicture(fileName, Range = Anchor)

		index_of_pic = app.ActiveDocument.InlineShapes.Count
		wdStory = 6

		app.Selection.HomeKey(Unit=wdStory)
		#app.ActiveDocument.InlineShapes(index_of_pic).Width = 28.35 * 5.5 
		#app.ActiveDocument.InlineShapes(index_of_pic).Height = 28.35 * 6.5

for title in output_hash:
	path ='c:\\Template-profile.doc'
	document = word.Documents.Open(  path  )

	word.Selection.Find.Execute('%raw_content%', False, False, False, False, False, True, 1, True, output_hash[ title  ]['raw_content'] , 2)
	word.Selection.Find.Execute('%raw_number%', False, False, False, False, False, True, 1, True, output_hash[ title  ]['raw_number'], 2)


	word.Selection.Find.Execute('%filterreads_content%', False, False, False, False, False, True, 1, True, output_hash[ title  ]['filterreads_content'], 2)
	word.Selection.Find.Execute('%filterreads_number%', False, False, False, False, False, True, 1, True, output_hash[ title  ]['filterreads_number'], 2)
	word.Selection.Find.Execute('%title%', False, False, False, False, False, True, 1, True, title, 2)
	word.Selection.Find.Execute('%Read1_title%', False, False, False, False, False, True, 1, True, title, 2)

	for each_staus in output_hash[ title  ]:
		if each_staus == 'raw':

			for each_tag in output_hash[ title  ][ each_staus  ]:
				if each_tag =='qualstats':
					MS_Word_Find_Past_Image(word, '#Read1_RQ#', output_hash[ title  ][ each_staus  ][ each_tag ])
				else:
					MS_Word_Find_Past_Image(word, '#Read1_RD#',output_hash[ title  ][ each_staus  ][ each_tag ])

		elif each_staus == 'filter':

			for each_tag in output_hash[ title  ][ each_staus  ]:
				if each_tag =='qualstats':
					MS_Word_Find_Past_Image(word, '#Read1_FQ#', output_hash[ title  ][ each_staus  ] [ each_tag ])
				else:
					MS_Word_Find_Past_Image(word, '#Read1_FD#', output_hash[ title  ][ each_staus  ] [ each_tag ])

	document.SaveAs(dir_output+'%s.doc'%(  title  ))
	document.Close()
