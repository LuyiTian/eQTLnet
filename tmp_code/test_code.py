#test_code.py

from motif import read_motif_PWM

motif_dict = {}
f_path = '/Users/luyi/data/motif/motifs.txt'
dist = []
for name,PWM in read_motif_PWM(f_path):
    motif_dict[name] = PWM


import pylab as pl 
pl.hist(dist, bins = 80)
pl.show()

'''
f_path = '/Users/luyi/data/motif/matches.txt'

num = 0
total_len = 0.
for line in open(f_path):
    items = line.strip().split()
    if items[1] == 'chr1':
        num +=1
        total_len += float(items[3])-float(items[2])
print num
print total_len
'''
