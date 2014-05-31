from math import sqrt,log
f = open("GD462.GeneQuantRPKM.50FN.samplename.resk10.txt",'r')
f_filter = open('EUR373.gene.cis.FDR5.best.rs137.txt','r')
f_out = open("test.txt",'w')
gene_out = open("gene.txt",'w')
sample_out = open("sample_id.txt",'w')
header = f.readline().strip().split('\t')
sample_out.write('\n'.join(header[4:]))
sample_out.close()
filter_list = []
for line in f_filter:
    filter_list.append(line.split('\t')[2])
f_filter.close()

def the_log(val):
    if val>0:return log(val)
    else:return -99999
p = 0
for line in f:
    all_vals = line.strip().split('\t')
    if all_vals[0] not in filter_list:continue
    tmp = ' '.join([str(the_log(float(val))) for val in all_vals[4:]])
    f_out.write(tmp+'\n')
    gene_out.write(str(p)+'\t'+'\t'.join(all_vals[:4])+'\n')
    p+=1
f.close()
f_out.close()
gene_out.close()

'''
import pylab as pl 
pl.hist(l, bins=50) 
pl.show()
''' 
