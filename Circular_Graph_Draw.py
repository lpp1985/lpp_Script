#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/11/12
"""

import PIL
from lpp import *
import subprocess
from PIL import Image,ImageDraw
from  jinja2 import FileSystemLoader,Environment
from Bio import SeqIO
from optparse import OptionParser
class Gene_Locus(object):
	#build block
	def __init__(self,start,end,cog="Other",strand=1):
		self.start = start
		self.stop = end
		self.cog = cog
		self.strand = strand
	def set_cog(cog):
		self.cog = cog
		
def circle_new(gc_view)  :
	ima = Image.open(gc_view).convert("RGBA")
	size = ima.size
	r2 = min(size[0], size[1])
	if size[0] != size[1]:
		ima = ima.resize((r2, r2), Image.ANTIALIAS)
	circle = Image.new('L', (r2, r2), 0)

	draw = ImageDraw.Draw(circle)
	draw.ellipse((500, 500, r2-500, r2-500), fill=255)
	alpha = Image.new('L', (r2, r2), 255)
	
	alpha.paste(circle, (0, 0))
	ima.putalpha(alpha)
	ima.save('test_circle.png')
	return ima
def ring_new(gene_view):
	ima = Image.open(gene_view).convert("RGBA")
	size = ima.size
	r2 = min(size[0], size[1])
	if size[0] != size[1]:
		ima = ima.resize((r2, r2), Image.ANTIALIAS)
	circle = Image.new('L', (r2, r2), 255)

	draw = ImageDraw.Draw(circle)
	draw.ellipse((5000, 4500, 15000, 15000), fill=0)
	alpha = Image.new('L', (r2, r2), 0)
	
	alpha.paste(circle, (0, 0))
	ima.putalpha(alpha)
	ima.save('test2_circle.png')
	return ima    
		



gview_root="/pub/SOFTWARE/Other/gview/"
cog_color = {"A":"rgb(255,0,0)",
             "B":"rgb(255,99,71)",
             "J":"rgb(240,128,128)",
             "K":"rgb(255,140,0)",
             "L":"rgb(255,20,147)",
             "D":"rgb(240,230,140)",
             "O":"rgb(189,183,107)",
             "M":"rgb(107,142,35)",
             "N":"rgb(34,139,34)",
             "P":"rgb(154,205,50)",
             "T":"rgb(50,205,50)",
             "U":"rgb(173,255,47)",
             "V":"rgb(0,250,154)",
             "W":"rgb(143,188,143)" ,
             "Y":"rgb(60,179,113)" ,
             "Z":"rgb(255,255,0)" ,
             "C":"rgb(0,255,255)" ,
             "G":"rgb(0,206,209)" ,
             "E":"rgb(70,130,180)" ,
             "F":"rgb(0,191,255)" ,
             "H":"rgb(0,0,255)" ,
             "I":"rgb(106,90,205)" ,
             "Q":"rgb(0,0,128)" ,
             "R":"rgb(190,190,190)" ,
             "S":"rgb(105,105,105)" ,
             "Unknown":"rgb(0,0,0)" ,
             "CDS":"rgb(0,0,153)" ,
             "tRNA":"rgb(153,0,0)" ,
             "rRNA":"rgb(153,0,153)" ,
             "Other":"rgb(51,51,51)" }




if __name__ == '__main__':
	usage = '''usage: python2.7 %prog [options] 
'''
	parser = OptionParser(usage =usage )



	parser.add_option("-g", "--GBK", action="store",
                              dest="gbk",
                              
                              help="GBK File")
	parser.add_option("-o", "--Output", action="store",
                              dest="output",
 
                              help="output png File")	
	(options, args) = parser.parse_args()
	gbk = options.gbk
	output = options.output
	GBK = SeqIO.parse( gbk,'genbank' )
	forward_gene = {}
	rever_gene = {}
	trna_gene = {}
	rrna_gene = {}
	for each_data in GBK:
		Name = each_data.name
		Length = "%s"%(len(each_data.seq))
		for each_feature in each_data.features[1:]:
			
			start,end = each_feature.location.start.real,each_feature.location.end.real
			#add cog on here
			if each_feature.type =="rRNA":
				data = Gene_Locus(start,end,strand=each_feature.strand,cog="rRNA")
				rrna_gene[data] = ''
				
			elif each_feature.type =="tRNA":
				data = Gene_Locus(start,end,strand=each_feature.strand,cog="tRNA")
				trna_gene[data] = ''
			############################
			elif "product" in each_feature.qualifiers:
				product = each_feature.qualifiers['product'][0]
				
				cate = re.findall("Category\[(\w)\]",product)
				if cate:
					cog = cate[0]
					if cog not in cog_color:
						print(cog)
					data = Gene_Locus(start,end,cog=cog)
			else:
				data = Gene_Locus(start,end)
			############################
			if each_feature.strand==1:
				forward_gene[data] = ''
			else:
				rever_gene[data] = ''
	
	
	
	
	
	templeloader = FileSystemLoader(gview_root)
	env = Environment(loader = templeloader)
	template = env.get_template('gview_templat.xml')
	END = open( "ring.xml" ,'w' )
	END.write(
		template.render(
	        {
	            "Accession":Name,
	            "Length":Length,
	            "DirCOG":forward_gene,
	            "RevCOG":rever_gene,
	            "cog_color":cog_color,
	            "rrna_gene":rrna_gene,
	            "trna_gene":trna_gene
	        }
		)
	)	
	END.close()
	command = "java -jar %s/gview.jar -t 100 -i %s -o ring.png -l circular -W 20000 -H 20000 -q high -f png"%(gview_root,END.name  )
	subprocess.check_output(command.split())
	command ="java -jar %s/gview.jar -t 100 -i %s -s %s -o circle.png -l circular  -W 12000 -H 12000 -q high -f png" %(gview_root,gbk,gview_root+'/myown.gss')
	subprocess.check_output(command.split())
	gc_graph = circle_new("circle.png")
	gene_view = ring_new("ring.png")

	new_img = Image.new("RGBA",gene_view.size, 255)
	
	

	new_img.paste(gene_view,mask=gene_view)
	new_img.paste(gc_graph,(4000,4000),mask=gc_graph)
	
	
	new_img.save(output)	