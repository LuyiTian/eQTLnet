import sqlite3
import cPickle
#from parameters import TF_peak_dir,data_dir,chrom_list
from read_processed_data import read_gene_mapping,get_gene_TSS_seg, read_expression_value
from overlap import cal_overlap
import os
class Db_SNP():
    def __init__(self, database_path):
        self.con = sqlite3.connect(database_path)
    def get_all_pos_by_chrm(self, chrm):
        cur=self.con.cursor()
        cur.execute("SELECT position FROM {chromosome}".format(chromosome = chrm))
        res = [i[0] for i in cur.fetchall()]
        print len(res)
        return res
    def get_pos_by_seg(self,seg,chrm):
        cur = self.con.cursor()
        cur.execute("SELECT * FROM {} WHERE position >= ? and position <= ?".format(chrm),seg)
        print len(cur.fetchall()[0])

if __name__ == "__main__":
    #gene_file_path = os.path.join('/Users/luyi/data',"Norm_raw_gene.txt")
    #id_to_gene, gene_to_id, _ = read_gene_mapping(gene_file_path)
    #for chrm in range(1,23):
    #    TSS_in_chr = get_gene_TSS_seg(gene_file_path,chrm)
    test_data = [(100000,500000),(1000000,1100000),(5000000,6000000),(6000000,7000000)]
    database_path = '/Users/luyi/data/Database/SNPdatabase.db'
    the_SNP = Db_SNP(database_path)
    the_SNP.get_pos_by_seg(test_data[0],'chr1')