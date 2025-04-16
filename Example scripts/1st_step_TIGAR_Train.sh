#!/bin/bash

chromnum=$chr_train


TIGARDIR=$path_to_tigar

DATADIR=$output_path

# TIGAR PARAMETERS
Gene_Exp_train_file=$path_to_expression_file
train_sample_path=$path_to_sampleID_file
geno_file_path=$path_to_genotype_file
genoformat=GT
mafval=0.01
hweval=0.00001
crossval=1

# OUT DIRECTORY
OUTDIR_DPR=${DATADIR}/DPR_Models
mkdir -p ${OUTDIR_DPR}

OUTDIR_EN=${DATADIR}/EN_Models
mkdir -p ${OUTDIR_EN}

# LOAD virtual environment for TIGAR
conda activate tigarenv

echo Running TIGAR train for chromosome "${chromnum}"


## Train DPR model
${TIGARDIR}/TIGAR_Model_Train.sh \
--model DPR \
--gene_exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_path} \
--genofile ${geno_file_path} \
--chr ${chromnum} \
--genofile_type vcf \
--format ${genoformat} \
--maf ${mafval} \
--hwe ${hweval} \
--cvR2 ${crossval} \
--dpr 1 \
--ES fixed \
--TIGAR_dir ${TIGARDIR} \
--out_dir ${OUTDIR_DPR} \
--thread ${SLURM_CPUS_PER_TASK}


## Train EN model
${TIGARDIR}/TIGAR_Model_Train.sh \
--model elastic_net \
--gene_exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_path} \
--genofile ${geno_file_path} \
--chr ${chromnum} \
--genofile_type vcf \
--format ${genoformat} \
--maf ${mafval} \
--hwe ${hweval} \
--cvR2 ${crossval} \
--TIGAR_dir ${TIGARDIR} \
--out_dir ${OUTDIR_EN} \
--thread ${SLURM_CPUS_PER_TASK}


# deactive virtual em
conda deactivate

