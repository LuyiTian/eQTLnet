import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from math import log,sqrt
x=[]
i = 0
for line in open("E:\\c_project\\read_from_file\\Debug\\"+'out_coexpression_all'):
    if i>=0:
        r = float(line.strip().split('\t')[2])
        if -1<r<1:
            x.append(0.5*log((1+r)/(1-r))*sqrt(420-3))
    i+=1
mu = 0
sigma = 1
num_bins = 50
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.5)
# add a 'best fit' line
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.xlabel('Fisher transformation')
plt.ylabel('Probability')
plt.title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()