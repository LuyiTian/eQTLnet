from parameters import processed_data_dir,chrom_list,snp_list_dir
from random import choice
from math import log
f = open(processed_data_dir+"direct_network",'r')
node_pair = []
p_values = []
p_MA = []
p_MB = []
p_BA = []
x=[]
y=[]
def cal_pro(MA, AB, M_B):
    return (1-MA)*(1-AB)*(M_B)
for line in f:
    items = line.strip().split('\t')
    edge = set([items[2], items[3]])
    x.append(float(items[6]))
    y.append(float(items[8]))
    if edge in node_pair:
        #x.append(cal_pro(float(items[5]), float(items[7]), float(items[8])))
        #y.append(cal_pro(p_MA[node_pair.index(edge)], p_BA[node_pair.index(edge)], p_values[node_pair.index(edge)]))
        pass
    else :
        node_pair.append(edge)
        p_MA.append(float(items[5]))
        p_MB.append(float(items[6]))
        p_BA.append(float(items[7]))
        p_values.append(float(items[8]))
z = [xi/yi for xi, yi in zip(x,y) if xi/yi>1.4 or xi/yi<0.6]
#new_x = [choice(z) for i in range(20000)]
#new_y = [choice(z) for i in range(20000)]


import pylab as pl 
pl.plot(x,y, 'x') 
pl.show()