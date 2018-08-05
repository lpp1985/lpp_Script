#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import re,sys,os
RAW = open(sys.argv[1],'rU')
path = os.path.dirname(os.path.abspath(sys.argv[1]) )+'/'
data_hash = {}
title = RAW.next().strip().split("\t")[1:]
for each_samp in title:
	data_hash[ each_samp ] = []
	
all_data = []

for line in RAW:
	line_l = line.split("\t")[1:]
	i=0
	for each_data in line_l:
		data_hash[ title[i] ].append(float(each_data))
		all_data.append(float(each_data))
		i+=1

	
data = [ data_hash[i]  for i in sorted(      data_hash,key = lambda x: int( re.search( "(\d+)\%$",x).group(1)  )      )]

sample = re.search("(^\S+)\-\d+\%$",i).group(1)
fig, ax1 = plt.subplots(figsize=(20,10),dpi=300)

fig.canvas.set_window_title('%s Rarefaction test'%(sample))
plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
bp = plt.boxplot(data,0,'')
plot_x = []
plot_y = []
for each_median in bp["medians"]:
	x,y = each_median.get_data()
	number = np.median( range(1,len('%.2f'%(y[0]) ) +1) )*-5 
	plt.annotate('%.2f'%(y[0]), xy=(np.average(x),y[0]),  xycoords='data',
		            xytext=(number, 20), textcoords='offset points',
		            arrowprops=dict(arrowstyle="->")
		            )	
	#plt.text(x[0],y[0]+5,"%.2f"%(y[0]))
	plot_x.append(np.average(x) )
	plot_y.append(np.average(y))
max_box = []
for each_box in bp["boxes"]:
	max_box.append(each_box.get_data()[1][2])
max_box = max(max_box)
#bp = plt.boxplot(data, notch=0, sym='-', vert=1, whis=1.5)
#plt.setp(bp['boxes'], color='black')
#plt.setp(bp['whiskers'], color='black')
#plt.setp(bp['fliers'], color='red', marker='+')

## Add a horizontal grid to the plot, but make it very light in color
## so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)
#ax1.set_yscale("log")
max_median = max([np.median(x)   for x in  data ])
#Hide these grid behind plot objects
ax1.set_axisbelow(True)
ax1.set_title('Rarefaction Testing')
ax1.set_xlabel('Data Perc %')
ax1.set_ylabel('RPKM')
plt.plot(plot_x, plot_y, 'g*')
plt.plot(plot_x, plot_y, 'r--')
xtickNames = plt.setp(ax1, xticklabels=[ "sample %s"%(i) for i in  sorted(data_hash  ,key = lambda x: int( re.search( "(\d+)\%$",x).group(1)  )) ])
ax1.set_ylim(  (0,max_box+20)  )  

#plt.show()
plt.savefig(path+"Rarefaction.png")
#plt.show()
