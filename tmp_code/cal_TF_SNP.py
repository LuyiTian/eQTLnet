'''
found SNPs in TF binding sites

TF_peak data is from ENCODE project:

'''
from parameters import TF_peak_dir,data_dir,chrom_list
from read_processed_data import read_gene_mapping,get_gene_TSS_seg, read_expression_value
from read_TF import read_narrowpeak
from overlap import cal_overlap
import os
#from numba import jit
from scipy.stats import linregress
def read_id_convert(id_convert_file_path):
    res_dict = {}
    id_convert_file = open(id_convert_file_path,'r')
    id_convert_file.readline()
    for line in id_convert_file:
        id1,id2 = line.strip().split('\t')
        #print id1,id2
        res_dict[id2] = id1
    return res_dict

def iter_snp(file_path):
    for line in open(file_path,'r'):
        yield line.strip().split(' ')

#@jit
def main():
    Thr = 0.05
    p_Thr = 0.0001
    snp_dir = os.path.join(data_dir,'snp_processed')
    gene_exp = read_expression_value(os.path.join(snp_dir,'tmp_exp_1'))#gene exp for each samples stored in gene_exp
    gene_file_path = os.path.join(data_dir,"Data/Norm_raw_gene.txt")
    id_to_gene, gene_to_id, _ = read_gene_mapping(gene_file_path)
    id_mapping = read_id_convert(os.path.join(data_dir,'Data/mart_export.txt'))
    out_file = open(os.path.join(data_dir,"Data/TF_SNP_res"),'w')
    cell_type = "Gm12878"
    TFdata_dir = os.path.join(TF_peak_dir, cell_type)
    file_list = os.listdir(TFdata_dir)
    TF_list = []
    for i in range(1,len(file_list)):
        TF = file_list[i].split("__")[1]
        print TF
        if TF in TF_list:continue
        TF_list.append(TF)
        if not id_mapping.has_key(TF):continue
        if not gene_to_id.has_key(id_mapping[TF]):continue
        TF_gene_id = gene_to_id[id_mapping[TF]]#eg: CTCF->2341
        all_TF_peak = read_narrowpeak(os.path.join(TFdata_dir,file_list[i]))
        for chrm in chrom_list:
            TSS_in_chr = get_gene_TSS_seg(gene_file_path,chrm)
            TF_peak_in_chr = [(it[1],it[2],ind) for ind,it in all_TF_peak.items() if it[0] == chrm]
            TF_peak_in_chr = sorted(TF_peak_in_chr, key = lambda x:x[0])
            if len(TF_peak_in_chr) == 0:
                print "cannot find TF in chromosome:",chrm
                continue
            tmp_dict = cal_overlap(TF_peak_in_chr,TSS_in_chr)
            if len(tmp_dict) == 0:
                print 'TF around genes not found in chromosome:',chrm
                continue
            tmp_list = []
            for vals in tmp_dict.values():
                tmp_list.extend([val for val in vals])
            tmp_list = tuple(set(tmp_list))
            TF_peak_in_chr = [it for it in TF_peak_in_chr if it[2] in tmp_list]
            TF_peak_in_chr = sorted(TF_peak_in_chr, key = lambda x:x[0])
            generator_snp = iter_snp(os.path.join(snp_dir,'tmp_snp_'+str(chrm)))
            snp_dict,TF_snp_dict = cal_overlap(generator_snp,TF_peak_in_chr)#SNP info for each samples stored in snp_dict 
            print len(snp_dict)
            gene_snp = {}
            for key, vals in tmp_dict.items():
                for val in vals:
                    if TF_snp_dict.has_key(val):
                        gene_snp.setdefault(key,[]).extend([s for s in TF_snp_dict[val]])
            for gene,snps in gene_snp.items():
                for snp in snps:
                    p_vals = [0.,0.,0.]
                    for typ in range(3):
                        if snp_dict[snp].count(typ)<len(snp_dict[snp])*Thr:
                            p_vals[typ] = 1.
                        else:
                            _, _, _, p_vals[typ], _ = linregress([the_exp for the_i,the_exp in enumerate(gene_exp[gene]) if snp_dict[snp][the_i] == typ],\
                                [the_exp for the_i,the_exp in enumerate(gene_exp[TF_gene_id]) if snp_dict[snp][the_i] == typ])
                    if min(p_vals)<= p_Thr:
                        pass
                        out_file.write(id_to_gene[gene]+'\t'+id_to_gene[TF_gene_id]+'\t'+str(chrm)+'\t'+str(snp)+'\t'+'\t'.join([str(the_p) for the_p in p_vals])+'\n')
    out_file.close()
if __name__ == "__main__":main()