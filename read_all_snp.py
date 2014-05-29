import multiprocessing
from io_snp import read_snp_by_chr
chrom_list = [13,14,15,16,17,18,19]
if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=4)
    result = []
    for i in chrom_list:
        msg = 'chr'+str(i)
        result.append(pool.apply_async(read_snp_by_chr, (msg, )))
    pool.close()
    pool.join()
    for res in result:
        print res.get()
    print "Sub-process(es) done."