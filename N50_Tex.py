#!/usr/bin/python
#coding:utf-8
# Author:   --<>
# Purpose: 
# Created: 2013/3/2
from lpp import *
RAW = open(sys.argv[1],'rU')
RAW.next()
data_l = RAW.next().strip().split("\t")
mean= data_l[5]
N50 = data_l[1]
TotalBase = data_l[-2]
BaseNumber = data_l[-1]
tex = r"""\begin{table}[H]
\begin{center}

    \caption{质控后Subreads统计表}
    
       \begin{tabularx}{\textwidth}{XcXX}
\hline
Mean Subread length&N50&Total Number of Bases&Number of Reads\\
 \hline
$\numprint{%s}&$\numprint{%s}&$\numprint{%s}&$\numprint{%s}\\
\hline
\end{tabularx}
\end{center}
\end{table}
"""%(mean,N50,TotalBase,BaseNumber )
END = open( sys.argv[2],'w' )
END.write(tex)
