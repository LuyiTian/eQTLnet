from copy import deepcopy
snp_list_dir = "E:\\snp_list\\"
def read_previous_eqtl(file_name):
    '''
    return dictionary:
    key: gene name i.e: ENSG00000137411.10
    value: position of the variant
    '''
    result_dict={}
    for line in open(file_name, 'r'):
        items = line.strip().split('\t')
        result_dict[items[2]] = float(items[6])#assume cis-qtl
        result_dict[float(items[6])] = items[2]
    return result_dict
def read_gene_mapping(file_name):
    '''
    id_to_gene:
        key: auto-increment id:0,1,2,3...
        value: gene name
    '''
    id_to_gene={}
    gene_to_id={}
    chr_to_gene={}
    for line in open(file_name,'r'):
        items = line.strip().split('\t')
        id_to_gene[int(items[0])] = items[1].split('.')[0]
        gene_to_id[items[1].split('.')[0]] = int(items[0])
        chr_to_gene.setdefault(items[3],[]).append(items[1].split('.')[0])
    return id_to_gene, gene_to_id, chr_to_gene    
def read_gene_sample_id(file_name):
    res_list = []
    for line in open(file_name):
        res_list.append(line.strip())
    return res_list
def read_undirected_n(file_name,id_to_gene, min_r = 0.0):
    res_dict={}
    for line in open(file_name,'r'):
        n1,n2 = [int(i) for i in line.strip().split('\t')[:2]]
        if -min_r <= float(line.strip().split('\t')[2]) <= min_r:continue
        res_dict.setdefault(id_to_gene[n1],[]).append(id_to_gene[n2])
        res_dict.setdefault(id_to_gene[n2],[]).append(id_to_gene[n1])
    return res_dict
def read_expression_value(file_name):
    res_list = []
    for line in open(file_name,'r'):
        res_list.append([float(i) for i in line.strip().split()])
    return res_list
def read_processed_snp(file_name):
    '''
    sample_id_list:ordered sample ids
    res_dict:
        key:position of variant 
            NOTE: positions like '89059044.5' exists, so use float() instead of int()
        value:a list of genotype, 0 for 0/0 1 for 1/0 2 for 1/1
    '''
    f = open(file_name,'r')
    sample_id_list = f.readline().strip().split('\t')
    res_dict={}
    for line in f:
        items = line.strip().split('\t')
        res_dict[float(items[0])] = [int(i) for i in items[1].split(',')]
    return sample_id, res_dict
def remap_expression(snp_sample_ids,gene_sample_ids,expression_matrix):
    '''
    expression list and snp list may have different sample orders
    force expression list to have the same order with snp list
    before: expression: s1, s3, s2. snp: s1, s2, s3
    after: expression: s1, s2, s3
    '''
    new_matrix=[[] for i in range(len(expression_matrix))]
    exclude_samples = []
    for the_id in snp_sample_ids:
        if the_id in gene_sample_ids:
            _i = gene_sample_ids.index(the_id)
            for ind,g in enumerate(expression_matrix):
                new_matrix[ind].append(g[_i])
        if the_id not in gene_sample_ids:
            print "WARNING: snp sample",the_id," is not in gene_sample_ids"
            exclude_samples.append(snp_sample_ids.index(the_id))
    return new_matrix, exclude_samples
def get_gene_TSS_seg(f_name,chr,distance = 3000):
    res = []
    for line in open(f_name,'r'):
        items = line.strip().split('\t')
        if items[3] == str(chr):
             res.append((int(items[4])-distance,int(items[4])+distance,int(items[0])))
    res = sorted(res,key = lambda x:x[0])
    return res
if __name__ == '__main__':
    file_name = "G:\\res_files\\EUR373.gene.cis.FDR5.best.rs137.txt"
    res = read_previous_eqtl(file_name)