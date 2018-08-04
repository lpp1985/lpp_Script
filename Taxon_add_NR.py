#!/usr/bin/env python

import os
import sys
from lpp import *
taxon_all = File_dict( open( sys.argv[1],'rU' )  ).read(  )

if __name__ == "__main__":
    unique_messages()
