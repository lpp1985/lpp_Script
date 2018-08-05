"""
Demo of the histogram (hist) function with a few features.

In addition to the basic histogram, this demo shows a few optional features:

    * Setting the number of data bins
    * The ``normed`` flag, which normalizes bin heights so that the integral of
      the histogram is 1. The resulting histogram is a probability density.
    * Setting the face color of the bars
    * Setting the opacity (alpha value).

"""
import sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

FILE = open( sys.argv[1],'rU'  )
data = []
for line in FILE:
	data.append( float(line) )
x = data

std = np.std(x)
mean = np.mean(x)
num_bins = 10
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
# add a 'best fit' line
y = mlab.normpdf(bins, mean, std)
plt.plot(bins, y, 'r--')
plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'Histogram of IQ: mean=%s$, srd=%s$'%( mean,std )
          )

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
pp = PdfPages(sys.argv[2])

plt.savefig(  pp ,format='pdf')
plt.show()
