from math import sqrt,log
from parameters import processed_data_dir, snp_list_dir
from copy import deepcopy
from scipy.stats import norm
from random import random as rd
def _log(val,offset = -0.001):
    if val<=offset:return -2.-rd()
    else:return log(val-offset)
def clean_log(val,MIN):
    if val <=0.:return MIN
    else:return log(val)
def read_exp(file_name,sample_list,gene_list = None,normalization = "Clean_Log"):
    f_out = open(processed_data_dir+normalization+"_raw_exp_value.txt",'w')
    gene_out = open(processed_data_dir+normalization+"_raw_gene.txt",'w')
    ##get sample from file
    sample_out = open(processed_data_dir+normalization+"_raw_sample_id.txt",'w')
    header = file_name.readline().strip().split('\t')
    sample_id = [i.split('.')[0] for i in header[4:] if i.split('.')[0] in sample_list]
    sample_index = [i for i,j in enumerate(header[4:]) if j.split('.')[0] in sample_list]
    sample_out.write('\n'.join(sample_id))
    sample_out.close()
    ##get gene and expression values in file
    tmp_dict ={}
    gene_index = 0
    norm_dist = []
    MIN = -99999.
    thr = 0.15
    for i in range(1,len(sample_index)+1):#get norm distribution cdf^-1(rank/(size+1))
        norm_dist.append(norm.ppf(float(i)/(len(sample_id)+1)))
    for line in file_name:
        all_vals = line.strip().split('\t')
        if (gene_list == None) or (all_val[1] in gene_list):
            gene_out.write(str(gene_index)+'\t'+'\t'.join(all_vals[:4])+'\n')
            tmp = [val for i, val in enumerate(all_vals[4:]) if i in sample_index]
            tmp_dict[gene_index] = [[i, val] for i,val in enumerate(tmp)]
        if normalization == 'Norm':
            tmp = deepcopy(tmp_dict[gene_index])
            tmp = sorted(tmp, key = lambda x:x[1])
            for i in range(len(tmp)):
                tmp[i][1] = norm_dist[i]
            for i, val in tmp:
                tmp_dict[gene_index][i][1] = val
        elif normalization == "Log":
            tmp_dict[gene_index] = [[i[0], _log(float(i[1]))] for i in tmp_dict[gene_index]]
        elif normalization == "Clean_Log":
            tmp_dict[gene_index] = [[i[0], clean_log(float(i[1]),MIN)] for i in tmp_dict[gene_index]]
            print len(tmp_dict[gene_index])
            if [i[1] for i in tmp_dict[gene_index]].count(MIN)>thr*len(tmp_dict[gene_index]):
                del tmp_dict[gene_index]
                continue
        else:
            f_out.write(' '.join([str(i[1]) for i in tmp_dict[gene_index]])+'\n')
            del tmp_dict[gene_index]
            gene_index+=1
        
    f_out.close()
    gene_out.close()
if __name__ == "__main__":
    f = open(processed_data_dir+"GD462.GeneQuantRPKM.50FN.samplename.resk10.txt",'r')
    f_snp =open(snp_list_dir+'res_chr1','r')
    all_sample = f_snp.readline()
    sample_list = all_sample.split('\t')
    print len(sample_list)
    f_snp.close()
    read_exp(f, sample_list,normalization = "Norm")