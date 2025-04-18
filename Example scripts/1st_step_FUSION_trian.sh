#!/bin/bash

chr=$chr_train

# Load modules 
module purge
module load R/4.2.2

# gene list
gene_list=$path_to_gene_list
n_row=`wc -l ${gene_list} | cut -d' ' -f 1`

for gene in $(seq 1 $n_row)
do

	line=`head -n $gene ${gene_list} | tail -n 1`
	gene_name=`echo ${line} | cut -d' ' -f 1`
  geno_dir=$path_to_genotype_per_gene
  
  cd ${geno_dir}
	b_dir=${gene_name} 
	tmp_dir=${gene_name}.tmp
	out_dir=${gene_name}
  cd $path_to_genotype_per_gene
  mkdir -p $path_to_genotype_per_gene/output

	Rscript /projects/YangLabData/qliu/TIGAR/input/FUSION/FUSION.compute_weights.R \
	--bfile $b_dir \
	--tmp $tmp_dir \
	--out $out_dir \
    --hsq_p 1.5 \
    --hsq_set 0.01\
	--PATH_plink $path_to_plink \
	--PATH_gcta $path_to_gcta \
	--PATH_gemma $path_to_gemma \
	--models top1,lasso,enet,blup \
    --save_hsq TRUE

done 