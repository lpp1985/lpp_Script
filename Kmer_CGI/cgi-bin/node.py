#!/usr/bin/python
import cgi
from pymongo import MongoClient
form = cgi.FieldStorage()


client = MongoClient('192.168.31.82', 27017)
node = form.getvalue('node')
org = form.getvalue('org')
db = client[ org ]
seq_db = db.Sequence
sequence = seq_db.find_one({"name":node})["seq"]
print "Content-type:text/html"
print
print """

<html>
<head>
        <title> Node  {node} in {org}</title>
<style>
table td{{
 white-space:pre-wrap;
 word-wrap: break-word;
}}
</style>
</head>        
       
<body>       
       <table cellpadding="0" width="100%" cellspacing="0" border="1"  style="word-break:break-all; word-wrap:break-all;">
<thead>
	<tr>
 		 <td bgcolor="#E3E2E2" align="center">
      			<font size=20 class="title3"><b>Sequence of Node <font  color="#FF0000">{node} </font> in<font  color="#5A95D"> {org}</font></b></font>
  		</td>
	</tr>
</thead>
<tbody>
	<tr>
	  	<td valign="bottom" align="left" width="100%"  style="table-layout:fixed">
  			 <br><font size=10>{sequence}</font></br>

  		</td>
	</tr>
</tbody>
</table>
</body>
</html>

""".format(node=node,org=org,sequence=sequence)
