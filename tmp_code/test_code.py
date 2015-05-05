#test_code.py

from motif import read_motif_PWM



motif_dict = {}
f_path = '/Users/luyi/data/motif/motifs.txt'
dist = []
for name,PWM in read_motif_PWM(f_path):
    motif_dict[name] = PWM
    dist.extend(PWM)


import pylab as pl 
pl.hist(dist, bins = 80)
pl.show()

'''
f_path = '/Users/luyi/data/motif/matches.txt'

res = []
for line in open(f_path):
    items = line.strip().split()
    if items[1] == 'chr1':
        res.append(items[0].split('_')[0])
res_list = [(the_id, res.count(the_id)) for the_id in list(set(res))]
res_list.sort(key=lambda x: x[1], reverse=True)

for i in res_list[:50]:
    print i
'''