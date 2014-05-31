snp_dir = "E:\\datasource\\"
#DEL_NUM = 42
import vcf
from parameters import snp_list_dir
def to_int(genotype):
    if genotype == "0|0" or genotype == "0/0": return "0"
    elif genotype == "1|0" or genotype =="1/0": return "1"
    elif genotype == "0|1" or genotype =="0/1": return "1"
    elif genotype == "1|1" or genotype =="1/1": return "2"
    elif genotype == None: return "-99999" 
    else: 
        print genotype
        raise IndexError
def read_snp_by_chr(chrom,p = 0.05):
    f = open(snp_list_dir+"res_"+chrom,'w')
    tmp_i = 0
    vcf_reader = vcf.Reader(open(snp_dir+"GEUVADIS"+'.'+chrom+'.'+'vcf.gz', "rb"))
    f.write('\t'.join(vcf_reader.samples)+'\n')
    for record in vcf_reader:
        tmp_l = [to_int(sample['GT']) for sample in record.samples]
        if tmp_l.count(0)<len(tmp_l)*(1-p):
            f.write(str(record.POS)+'\t'+','.join(tmp_l)+'\n')
            tmp_i+=1
        if tmp_i%10000 ==0: print "chrom:",chrom,tmp_i,"snps added"
    f.close()
    return "read "+chrom+" successfully"
if __name__ == "__main__":
    import sys
    chrom = 'chr'+sys.argv[1]
    read_snp_by_chr(chrom)
