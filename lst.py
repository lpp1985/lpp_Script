#!/usr/bin/env python
#coding:utf-8
# Author:  LPP --<lpp1985@hotmail.com>
# Purpose: 
# Created: 2011/9/19
from lpp import *
import os ,sys
class dir(object): 
	def __init__(self): 
		self.SPACE = "" 
		self.list = []
	
	def getCount(self, url):
		files = os.listdir(url)
		count = 0;
		for file in files:
			myfile = url + "\\" + file
			if os.path.isfile(myfile):
				count = count + 1
		return count
	def getDirList(self, url): 
		files = os.listdir(url) 
		fileNum = self.getCount(url)
		tmpNum = 0
		for file in files: 
			myfile = url + "/" + file 
			size = os.path.getsize(myfile) 
			if os.path.isfile(myfile): 
				tmpNum = tmpNum +1
				if (tmpNum != fileNum):
					self.list.append(str(self.SPACE) + "├─" + file + "\n")
				else:
					self.list.append(str(self.SPACE) + "└─" + file + "\n")
			if os.path.isdir(myfile): 
				self.list.append(str(self.SPACE) + "├─" + file + "\n") 
				# change into sub directory
				self.SPACE = self.SPACE + "│  " 
				self.getDirList(myfile) 
			    # if iterator of sub directory is finished, reduce "│  " 
				self.SPACE = self.SPACE[:-4] 
		return self.list 
	def writeList(self): 
		print(self.list) 
		print "ok" 

if __name__ == '__main__': 
	d = dir() 
	d.getDirList(sys.argv[1]) # input directory
	d.writeList() # write to file
