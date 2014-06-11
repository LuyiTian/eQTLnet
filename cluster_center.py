from Pycluster import kcluster
from parameters import processed_data_dir
from stat_coexp import get_in_out
def found_cluster_center(k = 100,repeat_times = 5):
    data = []
    for line in open(processed_data_dir+'Norm_raw_exp_value.txt'):
        vals = line.strip().split(' ')
        data.append([float(i) for i in vals])
    print len(data)
    for vals in data:
        if min(vals)<-8. or max(vals)>8.:
            print len(vals),min(vals),max(vals)
    belong,sum_val,nfound=kcluster(data,k,npass=repeat_times)
    print sum_val,nfound
    a_f = open(processed_data_dir+'gene_class.txt','w')
    for a,b in enumerate(belong):
        a_f.write(str(a)+'\t'+str(b)+'\n')
    a_f.close()
    res = get_in_out(processed_data_dir+'out_coexpression')
    final_select = [0.]*k
    tmp_biggest_len=[0]*k
    for ind,cla in enumerate(belong):
        if res.has_key(ind):
            if len(res[ind])>tmp_biggest_len[cla]:
                final_select[cla] = ind
                tmp_biggest_len[cla] = len(res[ind])
    for i in range(len(final_select)):
        if final_select[i] == 0.:
            print 'not found in class',i
            final_select[i] = random.choice(range(len(belong)))
    return final_select
def reduce_gene(final_select):
    gene_out = open(processed_data_dir+'subset_gene.txt','w')
    exp_out = open(processed_data_dir+'subset_exp_val.txt','w')
    gene_in = open(processed_data_dir+'Norm_raw_gene.txt','r')
    exp_in = open(processed_data_dir+'Norm_raw_exp_value.txt','r')
    ind = 0
    for ith, line in enumerate(gene_in):
        if ith in final_select:
            gene_out.write(str(ind)+'\t'+'\t'.join(line.strip().split('\t')[1:])+'\n')
            ind+=1
    for ith, line in enumerate(exp_in):
        if ith in final_select:
            exp_out.write(line)
if __name__ == "__main__":
    final_select = found_cluster_center()
    reduce_gene(final_select)
            