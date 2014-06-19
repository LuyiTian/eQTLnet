# distutils: language = c++
# distutils: sources = eqtl_m.cpp
cdef extern from "eqtl_m.h":
    cdef cppclass eQTL_Matrix:
        int smp_size
        int gene_num
        int SNP_size
        double *gene_exp
        double *SNP_chunk
        double *r_matrix
        eQTL_Matrix(int, int, int) except +
        void add_exp(double, int)
        void add_SNP(double, int)
        void reset_r()
        void exp_norm()
        void SNP_norm()
        void matrix_multi()

cdef class PyeQTL_Matrix:
    cdef eQTL_Matrix *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, int sample_size, int gene_number, int SNP_chunk_size):
        self.thisptr = new eQTL_Matrix(sample_size, gene_number, SNP_chunk_size)
    def __dealloc__(self):
        del self.thisptr
    def add_exp(self, double val, int i):
        self.thisptr.add_exp(val, i)
    def add_SNP(self, double snp, int i):
        self.thisptr.add_SNP(snp, i)
    def reset_r(self):
        self.thisptr.reset_r()
    def exp_norm(self):
        self.thisptr.exp_norm()
    def SNP_norm(self):
        self.thisptr.SNP_norm()
    def matrix_multi(self):
        self.thisptr.matrix_multi()
    def get_r(self,int i):
        return self.thisptr.r_matrix[i]

def eqtl_mapping(exp_file, snp_files, out_file, double thr_r, int snp_chunk_size):
    exp_res = []
    for line in open(exp_file,'r'):
        exp_res.append([float(i) for i in line.strip().split(' ')])
    sample_num = len(exp_res[0])
    gene_num = len(exp_res)
    print "sample number:",sample_num,'gene number:',gene_num
    M = PyeQTL_Matrix(sample_num, gene_num, snp_chunk_size)
    cdef int k
    cdef int j
    for k in range(gene_num):
        for j in range(sample_num):
            M.add_exp(exp_res[k][j],k*sample_num+j)
    M.exp_norm()
    del exp_res
    print 'successfully add and normalize expression value'
    the_out = open(out_file,'w')
    cdef double r
    for snp_file in snp_files:
        snp_f = open(snp_file,'r')
        print 'mapping snp file:',snp_f
        chunk = 0
        snp_name = ["."]*snp_chunk_size
        line = snp_f.readline()
        the_num = 0
        while line:
            items = line.strip().split(' ')
            snp_name[chunk] = items[0]
            #cdef int ind
            for ind in range(sample_num):
                M.add_SNP(float(items[ind+1]),chunk+ind*snp_chunk_size)
            if chunk == snp_chunk_size-1:
                M.SNP_norm()
                M.matrix_multi()
                #cdef int gene_i,snp_i
                for gene_i in range(gene_num):
                    for snp_i in range(snp_chunk_size):
                        r = M.get_r(gene_i*snp_chunk_size+snp_i)
                        if r>thr_r or r<-thr_r:
                            the_out.write(str(gene_i)+'\t'+snp_name[snp_i]+'\t'+str(r)+'\n')
                the_num+=1
                print the_num*snp_chunk_size,'snp processed'
                chunk = 0#reset chunk
                M.reset_r()#reset r_matrix
            chunk+=1
            line = snp_f.readline()
        snp_f.close()
    the_out.close()
         