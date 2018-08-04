#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2017/11/25
"""

from lpp import * 
def Rcode_write(dataset,file_prefix, format='pdf', colNum=100):
	'''generate R script for visualization'''
	ROUT = open(file_prefix + '.r','w')
	names=[]
	datas=[]
	for name, data in dataset:
		names.append(name)

		print >>ROUT, name + data

	tick_pos = [1,10,20,30,40,50,60,70,80,90,100]
	tick_lab = [1,10,20,30,40,50,60,70,80,90,100]

	# do not generate heatmap if only 1 sample
	if len(names) >=3:
		print >>ROUT, 'data_matrix' + ' <- matrix(c(' + ','.join(names) + '), byrow=T, ' +  'ncol=' + str(colNum) + ')'
		print >>ROUT, 'rowLabel <- c(' + ','.join(['"' + i + '"' for i in names]) + ')'
		print >>ROUT, '\n'
		print >>ROUT, '%s(\"%s.%s\")' % (format.lower(),file_prefix + ".heatMap",format.lower())
		print >>ROUT, 'rc <- cm.colors(ncol(data_matrix))'
		print >>ROUT, 'heatmap(data_matrix' + ', scale=c(\"none\"),keep.dendro=F, labRow = rowLabel ' + ',Colv = NA,Rowv = NA,labCol=NA,col=cm.colors(256),margins = c(6, 8),ColSideColors = rc,cexRow=1,cexCol=1,xlab="Gene body percentile (5\'->3\')", add.expr=x_axis_expr <- axis(side=1,at=c(%s),labels=c(%s)))' % (','.join([str(i) for i in tick_pos]), ','.join(['"' + str(i) + '"' for i in tick_lab]))
		print >>ROUT, 'dev.off()'


	print >>ROUT, '\n'

	print >>ROUT, '%s(\"%s.%s\")' % (format.lower(),file_prefix + ".curves",format.lower())	
	print >>ROUT, "x=1:100" 
	print >>ROUT, 'icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(%d)' % (len(names))
	if len(names) == 1:
		print >>ROUT, "plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0])

	elif  len(names) >=2 and len(names) <=6:  
		print >>ROUT, "plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0])
		for i in range(1,len(names)):
			print >>ROUT, "lines(x,%s,type='l',col=icolor[%d])" % (names[i], i+1)
		print >>ROUT, "legend(0,1,fill=icolor[%d:%d], legend=c(%s))" % (1,len(names), ','.join([ "'" + str(n) + "'" for n in names]))
	elif len(names) > 6:
		print >>ROUT, 'layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2), 4, 4, byrow = TRUE))'
		print >>ROUT, "plot(x,%s,type='l',xlab=\"Gene body percentile (5\'->3\')\", ylab=\"Coverage\",lwd=0.8,col=icolor[1])" % (names[0])
		for i in range(1,len(names)):
			print >>ROUT, "lines(x,%s,type='l',col=icolor[%d])" % (names[i], i+1)
		print >>ROUT, 'par(mar=c(1,0,2,1))'
		print >>ROUT, 'plot.new()'
		print >>ROUT, "legend(0,1,fill=icolor[%d:%d],legend=c(%s))" % (1,len(names), ','.join([ "'" + str(n) + "'" for n in names]))

	print >>ROUT, 'dev.off()'
	ROUT.close()
	os.system("Rscript %s" % (ROUT.name))



if __name__ == '__main__':
	all_data = sys.argv[1:-1]
	dataset = []
	for r_file in all_data:
		RAW = open (r_file, 'rU')
		line = RAW.next()
		dataset.append( line.split(" ", 1))
	Rcode_write(dataset, sys.argv[-1] )
