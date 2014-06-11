import os
import sys
from parameters import processed_data_dir, snp_list_dir, chrom_list
from read_processed_data import remap_expression,read_expression_value
import multiprocessing
#print sys.path[0]
#tmp = os.popen(sys.path[0]+'\\bin\\regression_clean D:\\gene.txt D:\\this.txt 4 3 0.1').readlines()
#print tmp
#for item in tmp:
#    print item.strip().split('\t')
def sample_match(exp_sample,snp_sample):
    index_list = []
    for i,sp_id in enumerate(snp_sample):
        if sp_id in exp_sample:
            index_list.append(i)
    return index_list
def eqtl_mapping(chr,exp_file,exp_sample_file,thr_r = 0.2):
    res = {}
    
    #print "process chr",chr
    
    snp_file = open(snp_list_dir+"res_chr"+str(chr),'r')
    snp_sample = snp_file.readline().strip().split('\t')
    exp_sample = [i.strip() for i in open(exp_sample_file)]
    filtered_index = sample_match(exp_sample, snp_sample)
    tmp_dir = "tmp"
    tmp_dir = os.path.join(processed_data_dir, tmp_dir)
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    expression_matrix = read_expression_value(exp_file)
    reordered_exp, _ = remap_expression([snp_sample[ind] for ind in filtered_index], exp_sample, expression_matrix)
    #if not os.path.isfile(tmp_dir+"\\tmp_exp_"+str(chr)):
    if True:
        tmp_exp = open(tmp_dir+"\\tmp_exp_"+str(chr),'w')
        for vals in reordered_exp:
            tmp_exp.write(' '.join([str(val) for val in vals])+'\n')
        tmp_exp.close()
    #print "remapping expression in chr",chr
    if not os.path.isfile(tmp_dir+"\\tmp_snp_"+str(chr)):
        tmp_snp = open(tmp_dir+"\\tmp_snp_"+str(chr),'w')
        print snp_file
        print tmp_snp
        for line in snp_file:
            snp_position, vals = line.strip().split('\t')
            vals = [j for i,j in enumerate(vals.split(',')) if i in filtered_index]
            tmp_snp.write(snp_position+' '+' '.join(vals)+'\n')
        tmp_snp.close()
        snp_file.close()
    '''
    #print 'filter snp samples in chr',chr
    
    cmd = "\\bin\\regression_clean "
    args = []
    args.append(tmp_dir+"\\tmp_exp_"+str(chr))#add gene expression file
    args.append(tmp_dir+"\\tmp_snp_"+str(chr))#add snp file
    args.append(str(len(filtered_index)))#add sample number
    args.append(str(len(reordered_exp)))#add gene number
    args.append(str(thr_r))#add thr_r
    #print "eqtl mapping in chr",chr
    print cmd+' '.join(args)
    for line in os.popen(sys.path[0]+cmd+' '.join(args)):
        print line
        snp_position, gene_index, r = line.strip().split('\t')
        res.setdefault(gene_index,[]).append((str(chr),snp_position,r))
    #os.remove(tmp_dir+"\\tmp_exp_"+str(chr))
    #os.remove(tmp_dir+"\\tmp_snp_"+str(chr))
    '''
    return res
if __name__ =="__main__":
    pool = multiprocessing.Pool(processes=4)
    exp_file = processed_data_dir+'Norm_raw_exp_value.txt'
    exp_sample_file = processed_data_dir+"Norm_raw_sample_id.txt"
    all_res = {}
    result = []
    for i in chrom_list:
        result.append(pool.apply_async(eqtl_mapping, (i,exp_file,exp_sample_file )))
    pool.close()
    pool.join()
    '''
    for res in result:
        for key,vals in res.get().items():
            all_res.setdefault(key,[]).extend(vals)
    out_file=open(processed_data_dir+"result_eqtl",'w')
    for key,vals in all_res.items():
        for val in vals:
            out_file.write(str(key)+'\t'+'\t'.join(val)+'\n')
    '''