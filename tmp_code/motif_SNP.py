import sqlite3
import cPickle
from parameters import TF_peak_dir,data_dir,chrom_list
from read_processed_data import read_gene_mapping,get_gene_TSS_seg, read_expression_value
from overlap import cal_overlap
import os
class Db_SNP():
    def __init__(self, database_path):
        self.con = sqlite3.connect(database_path)
    def get_pos_by_chrm(self, chrm):
        cur=self.con.cursor()
        cur.execute("SELECT position FROM {chromosome}".format(chromosome = chrm))
        res = [i[0] for i in cur.fetchall()]
        print len(res)
        return res

if __name__ == "__main__":
    gene_file_path = os.path.join(data_dir,"Data/Norm_raw_gene.txt")
    id_to_gene, gene_to_id, _ = read_gene_mapping(gene_file_path)
    for chrm in chrom_list:
        TSS_in_chr = get_gene_TSS_seg(gene_file_path,chrm)

    database_path = '/Users/luyi/eQTL/Database/SNPdatabase.db'
    the_SNP = Db_SNP(database_path)
    pos = the_SNP.get_pos_by_chrm('chr1')
    pos.sort()
    print pos[1],pos[-1]