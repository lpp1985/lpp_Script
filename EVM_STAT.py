#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/8/22
"""

from lpp import * 



if __name__ == '__main__':
	all_data = glob.glob(sys.argv[1])
	i = 0
	AUGUSTUS = open("augustus.list", 'w')
	AUGUSTUS.write("Gene\n")
	RNA =  open("RNA_SEQ.list", 'w')
	RNA.write("Gene\n")
	PRO =  open("Protein.list", 'w')
	PRO.write("Gene\n")
	SNAP = open("SNAP.list", 'w')
	SNAP.write("Gene\n")
	GENEMARK = open("GenMark.list",'w')
	GENEMARK.write("Gene\n")
	for e_f in all_data:
		
		RAW = open(e_f)
		all_block = RAW.read().split("\n\n")
		
		for e_b in  all_block:
			i += 1
			name = str(i)
			if "SNAP" in e_b:
				SNAP.write(name + '\n')
			if "AUGUSTUS" in e_b:
				AUGUSTUS .write(name + '\n')
			if "GeneMark" in e_b:
				GENEMARK.write(name + '\n')
			if "assembler" in e_b:
				RNA.write(name + '\n')
			if "GeneWise" in e_b:
				PRO.write(name + '\n' )
			
	
	VENN_R = open( "Draw.R",'w'   )
	VENN_R.write("""
require("VennDiagram")
temp = venn.diagram(
	 x = list(

	    """)
		
		

	end_list = []
	file_hash = {"EST": RNA.name, "Protein": PRO.name, "SNAP": SNAP.name, "GeneMarks": GENEMARK.name, "AUGUSTUS": AUGUSTUS.name,}
	for category, file_name in file_hash.items():
		
		end_list.append(  """    %s =  read.delim( \"%s\", header=TRUE, stringsAsFactors=TRUE )$Gene"""%(
		    category,
		    file_name
		)
		                  )
	VENN_R.write(",".join(end_list))

	VENN_R.write("""),
        filename = NULL,
    
        fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        alpha = 0.50,
    
        cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 
        1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
        margin = 0.05,
        cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.cex = 1,
    
        ) 
    pdf("%s")
    grid.draw(temp)    
    dev.off()  
    tiff( "%s"  )  
    grid.draw(temp) 
    dev.off()
        """%(
               'stat.pdf',
               'stat.tiff'
    
    
           )
    
    
           )	
	
		
