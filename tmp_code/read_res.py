from parameters import data_dir
import os
import matplotlib.pyplot as plt
from read_processed_data import read_gene_mapping,get_gene_TSS_seg, read_expression_value
from cal_TF_SNP import iter_snp
'''
out_f = open(os.path.join(data_dir,"Data/reduce_TF_SNP_res"),'w')
p_lo = 0.000001
p_hi = 0.1
for line in open(os.path.join(data_dir,"Data/TF_SNP_res"),'r'):
    p_list = [float(i) for i in line.strip().split('\t')[4:]]
    for i in p_list:
        if p_hi<i<1. and min(p_list)<p_lo:
            out_f.write(line)
            continue
out_f.close()
'''


snp_dir = os.path.join(data_dir,'snp_processed')
gene_exp = read_expression_value(os.path.join(snp_dir,'tmp_exp_1'))#gene exp for each samples stored in gene_exp
gene_file_path = os.path.join(data_dir,"Data/Norm_raw_gene.txt")
_, gene_to_id, _ = read_gene_mapping(gene_file_path)
convert_id = {}
for key, val in gene_to_id.items():
    convert_id[key.split('.')[0]] = val


l = "ENSG00000188976    ENSG00000154727 1   896064.0    3.69173268209e-08   0.375896279147  1.0"


items = l.split()
gene1 = items[0]
id_1 = convert_id[gene1]
gene2 = items[1]
id_2 = convert_id[gene2]
chrm = items[2]
position = float(items[3])

generator_snp = iter_snp(os.path.join(snp_dir,'tmp_snp_'+str(chrm)))
for items in generator_snp:
    if float(items[0]) == position:
        snp_list = [float(it) for it in items[1:]]
        print 'find SNP'
        break

fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True,sharey=True)
ax = axs[0]
ax.plot([x for i,x in enumerate(gene_exp[id_1]) if snp_list[i] == 0 ], [y for i,y in enumerate(gene_exp[id_2]) if snp_list[i] == 0 ], 'x')
ax.set_title('0|0')

ax = axs[1]
ax.plot([x for i,x in enumerate(gene_exp[id_1]) if snp_list[i] == 1 ], [y for i,y in enumerate(gene_exp[id_2]) if snp_list[i] == 1 ], 'x')
ax.set_title('1|0')

ax = axs[2]
ax.plot([x for i,x in enumerate(gene_exp[id_1]) if snp_list[i] == 2 ], [y for i,y in enumerate(gene_exp[id_2]) if snp_list[i] == 2 ], 'x')
ax.set_title('1|1')
fig.suptitle("x: %s  y: %s  SNP: %d %d" % (gene1,gene2,int(chrm),int(position)))
plt.show()