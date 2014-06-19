class eQTL_Matrix
{
    public:
        int smp_size;
        int gene_num;
        int SNP_size;
        double *gene_exp;
        double *SNP_chunk;
        double *r_matrix;
        eQTL_Matrix(int sample_size, int gene_number, int SNP_chunk_size);
        ~eQTL_Matrix();
        void add_exp(double val,int i);
        void add_SNP(double snp,int i);
        void reset_r();
        void exp_norm();//normalize data so they have mean=0; sum of square = 1
        void SNP_norm();//normalize data so they have mean=0; sum of square = 1
        void matrix_multi();//TODO: find quicker way to do multiplication 
};