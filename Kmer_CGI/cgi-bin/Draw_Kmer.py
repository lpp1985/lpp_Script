#!/home/user/miniconda2/bin/python
import matplotlib,sys,re
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
import numpy as np 
from matplotlib.patches import Circle,Rectangle,Polygon
from matplotlib.patheffects import withStroke
from pymongo import MongoClient
from collections import defaultdict,namedtuple
client = MongoClient('192.168.31.82', 27017)
db = client[sys.argv[2]]

class fasta_check(object):
	def __init__(self,file_handle):
		assert isinstance( file_handle,file ),'The Input paramater must be a File Handle'
		self.file=iter(file_handle)
		for line in self.file:
			if line[0]=='>':
				self.define=line
				break
	def __iter__(self):
		return self
	def next(self):
		if not self.define:
			raise StopIteration

		name=self.define
		self.define=''        
		s=[]
		for line in self.file:
			if line[0]!='>':
				s.append(line)
			else:
				self.define=line
				break
		s=''.join(s)
		return (name,s)




class Ddict(defaultdict,dict):
	def __init__(self):
		defaultdict.__init__(self, Ddict)
	def __repr__(self):
		return dict.__repr__(self)
seq_db = db.Sequence
loc_db = db.Location
loc_db.delete_many({})
all_nodes ={}
graph = Ddict()
hierachy= {}
node_length = {}
each_dminion_max = {}
RAW = fasta_check( open(sys.argv[1],'rU'))
max_length = 0
for t,s in RAW:
	s = re.sub("\s+","",s)
	if len(s) > max_length:
		max_length=len(s)
	name_l = t[1:].strip().split("; ")

	for i in xrange(1,len(name_l)):
		graph[ name_l[i-1]  ][ name_l[i]  ]=""
		if name_l[i-1] not in all_nodes:
			all_nodes[ name_l[i-1]   ] = ""
			if i-1 not in hierachy:
				hierachy[i-1]=[name_l[i-1]   ]
			else:
				hierachy[i-1].append( name_l[i-1]   )
		if name_l[i] not in all_nodes:    
			if i not in hierachy:
				hierachy[i]=[name_l[i]   ]
			else:
				hierachy[i].append( name_l[i]   )            
			all_nodes[ name_l[i]   ] = ""

max_num = 0
for i in hierachy:
	if i not in each_dminion_max:
		each_dminion_max[i]=0
	if len(hierachy[i])>max_num:
		max_num=len(  hierachy[i] )
	for each_node in hierachy[i]:
		length = len(seq_db.find_one({"name":each_node})["seq"])
		node_length[ each_node ] = length
		if length >each_dminion_max[i]:
			each_dminion_max[i]=length





def Million(x, pos):

	return "%d K" % (x/1e3)

def Horizion_Line( start_l,end_l   ):
	ax.plot( start_l, end_l, linewidth=0,zorder=5,linestyle="-.",color='r')

def Vertical_Line( start_l,end_l   ):
	ax.plot( start_l, end_l, linewidth=0,zorder=5,linestyle="-",color='b')    
def TransCoords(x,y,coord="x"):
	inv = ax.transData.inverted()
	start,end = inv.transform( ax.transAxes.transform([x,y]) )

	if coord=="x":
		return start,y
	else:
		return x,end

def Text_Draw( m1 ,  m2 ,text, coord="x",color="red",roat=90):
	x1,x2 = m1
	y1,y2 = m2
	if coord=="x":
		y_loc =np.average((y1,y2))
		x_loc = TransCoords(x1, y2,"x")[0]
	elif coord=="y":
		x_loc =np.average((x1,x2))

		y_loc = TransCoords(x1, y2,"y")   [1]     
	else:
		y_loc =np.average((y1,y2))
		x_loc =np.average((x1,x2))

	ax.text(x_loc,y_loc,text,backgroundcolor="none",
            ha='center', va='center', weight='bold', zorder=1,color=color,clip_on=False,fontsize=10,rotation=roat)   

def Rectangle_Draw(text,x,y,width,height  ):

	sequence = seq_db.find_one({"name":text})["seq"]
	rec = Rectangle((x, y), width,height, clip_on=False, zorder=0, linewidth=1,
                    edgecolor='black', facecolor=(0.22, 0.41, 1, .6),
                    path_effects=[withStroke(linewidth=0.5, foreground='w')],url="/Kmer/node.py?node=%s&org=%s"%(text.replace("+","%2B"),sys.argv[2]))
	start =  ax.transAxes.transform((x,y))

	end = ax.transAxes.transform((x+width,y+height))
	ax.add_artist(rec)    
	LOC = {
            "name":"rec",
    "location":[start[0],start[1],end[0],end[1]]
    }
	loc_db.insert(LOC)
	return ((x,x+width,y,y+height),  (x,y+height/2),  (x+width,y+height/2))# Output Text position,left Posion,Right Position




##Draw Figure and Clearfy unit of Axe
fig = plt.figure(figsize=(40, 40),dpi=300)
ax = fig.add_subplot(1, 1, 1, aspect=1)


height = 8000

layer = max(hierachy)


## X axis Setting
ax.xaxis.set_major_locator(MultipleLocator(5e3))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_major_formatter( FuncFormatter(Million)  )
ax.set_xlim(0, max_length*3)

## X axis Setting
ax.yaxis.set_major_locator(MultipleLocator(5e6))
ax.yaxis.set_major_formatter( FuncFormatter(Million)  )
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

ax.set_ylim(0, height)


distance,rec_height  = ax.transData.inverted().transform( ax.transAxes.transform([0.01,0.02]) )

each_node_coord = {}


# Draw Each Layer            
e_start = 50
for each_layer in hierachy:
	step = height/(len(hierachy[each_layer])+1)
	each_layer_hor = step
	for node in hierachy[each_layer]:

		each_node_coord[node]=Rectangle_Draw(node, e_start,each_layer_hor ,node_length[node],rec_height, )
		x_start,x_end,y_start,y_end = each_node_coord[node][0]
		Text_Draw((x_start,x_end), (y_start,y_end), node,coord="no",color="black")
		each_layer_hor +=step
	e_start +=each_dminion_max[ each_layer ]+distance

##Draw Line
for start in graph:
	for end in graph[start]:
		lef_loc = each_node_coord[start ][-1]
		right_loc = each_node_coord[end ][1]


		line_poly = Polygon( [ [lef_loc[0],  lef_loc[1]],[lef_loc[0],  lef_loc[1]-20],[right_loc[0],  right_loc[1]-20],[right_loc[0],  right_loc[1]] ],linewidth=1,zorder=5,color="black",url="/Kmer/Path.py?start=%s&end=%s&org=%s"%(start.replace("+","%2B"), end.replace("+","%2B"),sys.argv[2]) )
		ax.add_artist(line_poly)


#


plt.savefig(sys.argv[3],dpi=300)
#plt.show()
#plt.show()
