#!/usr/bin/python
import cgi
from pymongo import MongoClient
form = cgi.FieldStorage()


client = MongoClient('192.168.31.82', 27017)
start = form.getvalue('start')
end = form.getvalue('end')
kmer =  form.getvalue('kmer')
if not kmer:
    kmer = 41
else:
    kmer=int(kmer)
org = form.getvalue('org')
db = client[ org]
seq_db = db.Sequence
sequence = """<font color="#00FF00">"""+seq_db.find_one({"name":start})["seq"]+"</font>"
sequence+="""<font color="#FF0000">N</font>"""
sequence += """<font color="#0000F0">"""+seq_db.find_one({"name":end})["seq"][kmer-1:]+"</font>"
print "Content-type:text/html"
print
print """

<html>
<head>
        <title> Node  {start} -->{end} in {org}</title>
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
      			<font size=20 class="title3"><b>Sequence of Node  <font color="#00FF00"> {start}</font> --> <font color="#0000F0">{end}</font> in<font  color="#5A95D"> {org}</font></b></font>
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

""".format(start=start,end=end,org=org,sequence=sequence)
