#!/usr/bin/python
# -*- coding: UTF-8 -*-
from pymongo import MongoClient
client = MongoClient('192.168.31.82', 27017)
db = client["C5796"]
loc_db = db.Location
print "Content-type:text/html"
print                               # 空行，告诉服务器结束头部
print '<html>'
print '<head>'
print '<meta charset="utf-8">'
print '<title>Hello World - 我的第一个 LPP CGI 程序！</title>'
print '</head>'
print '<body>'
print '<h2>Hello World! LPP 我是来自菜鸟教程的第一CGI程序</h2>'
print """<img src="/Graph/test.png" name="pathwayimage" usemap="#mapdata" border="1">"""

print("""<map name="mapdata">""")
for i in loc_db.find({'name': "rec"}):
    start_x,start_y,end_x,end_y = i["location"]
    print("""<area shape="rect" coords="%d,%d,%d,%d" href="/dbget-bin/www_bget?K05774+2.7.4.23+R06836" title="K05774 (phnN), 2.7.4.23, R06836">"""%(start_x,start_y,end_x,end_y ))
print("""<area shape="rect" coords="%d,%d,%d,%d" href="/dbget-bin/www_bget?K05774+2.7.4.23+R06836" title="K05774 (phnN), 2.7.4.23, R06836">"""%(0,0,20,20 ))
print("""</map>""")
print '</body>'
print '</html>'
