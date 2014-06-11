#
#using method described in Aten et al(2007)
#
#Reference:
#Aten JE, Fuller TF, Lusis AJ, Horvath S. 
#Using genetic markers to orient the edges in quantitative trait networks: the NEO software. BMC Syst Biol. 
#2008;2:34. doi:10.1186/1752-0509-2-34.
#
#
#####################################################
from read_processed_data import *
from scipy.stats import norm, linregress
from parameters import processed_data_dir,chrom_list,snp_list_dir
from math import sqrt,log
#norm.cdf(1.96)
#slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

def cal_partial_correlation_r(node_M, node_A, node_B, lens):
    _, _, r_MA, p_MA, _ = linregress(node_M,node_A)
    _, _, r_MB, p_MB, _ = linregress(node_M,node_B)
    _, _, r_BA, p_BA, _ = linregress(node_A,node_B)
    partial_r = (r_MB - r_MA*r_BA)/sqrt((1-r_MA**2)*(1-r_BA**2))
    Z_score = abs(0.5*sqrt(lens-3)*log((partial_r+1)/(1-partial_r)))#Fisher's Z transform
    p_value = 2*norm.cdf(-Z_score)
    return partial_r, p_value,p_MA, p_MB, p_BA
    
if __name__ == "__main__":
    #########
    eqtl_name = processed_data_dir+"EUR373.gene.cis.FDR5.best.rs137.txt"
    gene_expression_name = processed_data_dir+"gene_expression.txt"
    gene_name = processed_data_dir+"gene.txt"
    gene_sample_name = processed_data_dir+"sample_id.txt"
    undi_network_name = processed_data_dir+"out_coexpression"
    #########
    eqtl_gene_snp = read_previous_eqtl(eqtl_name)
    id_to_gene, gene_to_id, chr_to_gene = read_gene_mapping(gene_name)
    gene_sample_ids = read_gene_sample_id(gene_sample_name)
    network_dict = read_undirected_n(undi_network_name,id_to_gene)
    expression_matrix = read_expression_value(gene_expression_name)
    old_sample_id = []
    new_matrix = []
    result_file = open(processed_data_dir+"direct_network",'w')
    print "start"
    for chrom in chrom_list:
        genes_in_chr = chr_to_gene[str(chrom)]
        print "start chr:",chrom
        snp_f = open(snp_list_dir+"res_chr"+str(chrom),'r')#read all SNP in this chromosome
        sample_id = snp_f.readline().strip().split('\t')
        if old_sample_id == sample_id:
            pass
        else:
            #if sample orders are different from previous one, change orders
            old_sample_id = sample_id
            new_matrix, exclude_samples = remap_expression(sample_id,gene_sample_ids,expression_matrix)
        snp_in_chr = [eqtl_gene_snp[gene] for gene in genes_in_chr if eqtl_gene_snp.has_key(gene)]
        print len(genes_in_chr),'genes in chromosome.',len(snp_in_chr),'snps in chromosome'
        for line in snp_f:
            items = line.strip().split('\t')
            snp = float(items[0]) 
            genotypes = [int(locus) for i, locus in enumerate(items[1].split(',')) if i not in exclude_samples]
            if snp in snp_in_chr:
                the_gene = eqtl_gene_snp[snp]
                if not network_dict.has_key(the_gene):
                    #print 'no gene is associated to',the_gene
                    continue
                for associated_gene in network_dict[the_gene]:
                    M = []
                    A = []
                    B = []
                    for i in range(len(genotypes)):
                        if genotypes[i]>-999 and new_matrix[gene_to_id[the_gene]]>-999 and new_matrix[gene_to_id[associated_gene]]>-999:
                            M.append(genotypes[i])
                            A.append(new_matrix[gene_to_id[the_gene]][i])
                            B.append(new_matrix[gene_to_id[associated_gene]][i])
                    r, p_value, p_MA, p_MB, p_BA = cal_partial_correlation_r(M, A, B,len(M))
                    result_file.write(str(chrom)+'\t'+str(snp)+'\t'+the_gene+'\t'+associated_gene+\
                    '\t'+str(r)+'\t'+str(p_MA)+'\t'+str(p_MB)+'\t'+str(p_BA)+'\t'+str(p_value)+'\n')
        snp_f.close()
    result_file.close()