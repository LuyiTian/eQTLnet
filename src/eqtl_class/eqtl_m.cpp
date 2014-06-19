#include "eqtl_m.h"
#include <math.h>
eQTL_Matrix::eQTL_Matrix(int sample_size, int gene_number, int SNP_chunk_size)
{
    smp_size = sample_size;
    gene_num = gene_number;
    SNP_size = SNP_chunk_size;
    gene_exp = new double[sample_size*gene_number];
    SNP_chunk = new double[sample_size*SNP_chunk_size];
    r_matrix = new double[gene_number*SNP_chunk_size];
    for (int i=0;i<gene_number*SNP_chunk_size;i++)
    {
        r_matrix[i] = 0;
    }
}
eQTL_Matrix::~eQTL_Matrix()
{
    delete []gene_exp;
    delete []SNP_chunk;//here SNP chunk is the transpose of SNP_matrix. so row is sample, col is snp
    delete []r_matrix;
}
void eQTL_Matrix::add_exp(double val,int i)
{
    gene_exp[i] = val;
}
void eQTL_Matrix::add_SNP(double SNP,int i)
{
    SNP_chunk[i] = SNP;
}
void eQTL_Matrix::reset_r()
{
    for (int i=0;i<gene_num*SNP_size;i++)
    {
        r_matrix[i] = 0;
    }
}
void eQTL_Matrix::exp_norm()
{
    double mean ;
    double sum_of_s;
    for(int i=0;i<gene_num;i++)
    {
        mean = 0;
        sum_of_s = 0;
        for (int j=0;j<smp_size;j++)mean+=gene_exp[i*smp_size+j];
        mean=mean/smp_size;
        for (int j=0;j<smp_size;j++)
        {
            gene_exp[i*smp_size+j] = gene_exp[i*smp_size+j]-mean;
            sum_of_s+=gene_exp[i*smp_size+j]*gene_exp[i*smp_size+j];
        }
        sum_of_s = sqrt(sum_of_s);
        for (int j=0;j<smp_size;j++) gene_exp[i*smp_size+j] = gene_exp[i*smp_size+j]/sum_of_s;
    }
}
void eQTL_Matrix::SNP_norm()
{
    double mean ;
    double sum_of_s;
    for (int i=0;i<SNP_size;i++)
    {
        mean = 0;
        sum_of_s = 0;
        for (int j=0;j<smp_size;j++)mean+=SNP_chunk[i+j*SNP_size];
        mean=mean/smp_size;
        for (int j=0;j<smp_size;j++)
        {
            SNP_chunk[i+j*SNP_size] = SNP_chunk[i+j*SNP_size]-mean;
            sum_of_s+=SNP_chunk[i+j*SNP_size]*SNP_chunk[i+j*SNP_size];
        }
        sum_of_s = sqrt(sum_of_s);
        for (int j=0;j<smp_size;j++) SNP_chunk[i+j*SNP_size] = SNP_chunk[i+j*SNP_size]/sum_of_s;
    }
}
void eQTL_Matrix::matrix_multi()
{
    double tmp;
    for(int i=0;i<gene_num;i++)
    for(int k=0;k<smp_size;k++)
    {
        tmp=gene_exp[i*smp_size+k];
        for(int j=0;j<SNP_size;j++)
            r_matrix[i*SNP_size+j]+=tmp*SNP_chunk[k*SNP_size+j];
    }
}