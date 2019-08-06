#!/usr/bin/python
import pandas as pd
import sys

tb = pd.read_csv(sys.argv[1])
tb.to_csv( sys.argv[2],sep="\t",index=0  )
