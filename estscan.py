#!/usr/bin/env python
#coding:utf-8
"""
  Author:   --<>
  Purpose: 
  Created: 2014/9/15
"""

import urllib,urllib2,re, HTMLParser,time
from lpp import *
RAW = fasta_check(open(sys.argv[1],'rU'))
def get_data(data,result):
    url = "http://myhits.isb-sib.ch/cgi-bin/estscan"
    values = {
        "species":"Drosophila_melanogaster.smat",
        "text":data,
       "action":"ESTScan",
       "indelpenalty":-50
              }
    data = urllib.urlencode(values)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    Data =  response.read()
    all_need = re.findall("""<td class='desc_summary' valign='top'><input type='radio' name='seqname' value='([^\']+)'""",Data,re.MULTILINE)
    html_parse = HTMLParser.HTMLParser()
    END=open(result,'a')
    all_need = map(lambda x:html_parse.unescape(x), all_need)
    all_need = map(lambda x:re.sub('<[^>]+>','',x),all_need)
    
    #all_need = filter(lambda x:"' onclick='theForm.text.value=this.value;' />" in x, all_need)
    END.write('\n'.join(all_need) +'\n')
num=0
cache = ""
end = sys.argv[2]
for t,s in RAW:
    num+=1
    cache+=t+s
    if num%50==0:
        get_data(cache, end)
        cache = ""
        #time.sleep(1)
else:
    get_data(cache, end)
