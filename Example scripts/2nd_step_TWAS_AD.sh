#!/bin/bash


chr=$chr_test

## GWAS Zscore files
Zscore=$path_to_gwas_file
## LD covariance file
LD_file=$path_to_LD_file

# TIGAR DIRECTORY
TIGAR_dir=$path_to_tigar

# LOAD virtual environment for TIGAR
conda activate tigarenv

############# RUN TWAS

# TIGAR DPR weights
weight=$path_to_DPR_weight
gene_anno=$path_to_gene_annotation_file
out_dir=$path_to_DPR_output

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${gene_anno} \
--Zscore ${Zscore} \
--weight ${weight} \
--LD ${LD_file} \
--chr ${chr} \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir} \
--thread ${SLURM_CPUS_PER_TASK}

# TIGAR EN weights
weight=$path_to_EN_weight
gene_anno=$path_to_gene_annotation_file
out_dir=$path_to_EN_output

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${gene_anno} \
--Zscore ${Zscore} \
--weight ${weight} \
--LD ${LD_file} \
--chr ${chr} \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir} \
--thread ${SLURM_CPUS_PER_TASK}

# TIGAR FUSION weights
weight=$path_to_FUSION_weight
gene_anno=$path_to_gene_annotation_file
out_dir=$path_to_FUSION_output

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${gene_anno} \
--Zscore ${Zscore} \
--weight ${weight} \
--LD ${LD_file} \
--chr ${chr} \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir} \
--thread ${SLURM_CPUS_PER_TASK}

# deactive virtual em
conda deactivate

