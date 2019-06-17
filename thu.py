#!/usr/bin/python
from PIL import  Image
import sys
Image.MAX_IMAGE_PIXELS = 1000000000
infile = sys.argv[1]
outfile = sys.argv[1]+".thu.png"
im = Image.open(infile)

out = im.resize((4000,4000),Image.ANTIALIAS) #resize image with high-quality
out.save(outfile)
