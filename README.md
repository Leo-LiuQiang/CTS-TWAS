# CTS-TWAS: Cell-type-speciffic Transcriptomic-wide Association Studies
### An omnibus transcriptomic-wide association study (TWAS) using cell-type-specific snRNA-seq data and integrating multiple statistical methods

![image](https://github.com/Leo-LiuQiang/CTS-TWAS/blob/main/xWAS-O framework.png)

## How It's Made:

### Reference

# Environment Setup
### 1. Download [TIGAR](https://github.com/yanglab-emory/TIGAR) and complete its software setup
 
* [BGZIP](http://www.htslib.org/doc/bgzip.html)
* [TABIX](http://www.htslib.org/doc/tabix.html) 
* Python 3.5 modules/libraries: pandas, numpy, scipy, sklearn, statsmodels
	* creating environment using conda
	```
	# create the environment tigarenv
	conda create --name tigarenv python=3.5 pandas numpy scipy scikit-learn statsmodels
	# deactivate the conda environment
	conda deactivate
	# activate the environment
	conda activate tigarenv
	# set the PYTHONPATH
	export PYTHONPATH=${CONDA_PREFIX}/lib/python3.5/site-packages/:$PYTHONPATH
	```
	* creating environment using pip
	```
	# install pip
	# install virtualenv
	python3 -m pip install --user virtualenv
	# cd to preferred install_directory
	cd ${install_dir}
	# create the virtual environment tigarenv in the current directory
	python3 -m virtualenv tigarenv --python=python3.5
	# activate the environment
	source ${install_dir}/tigarenv/bin/activate
	# install the packages
	python3 -m pip install numpy==1.15.2 pandas==0.23.4 scikit-learn==0.20.0 scipy==1.1.0 statsmodels==0.9.0
	# deactivate the environment
	deactivate
 	```
### 2. Download [FUSION tool](http://gusevlab.org/projects/fusion/) and complete its software setup

* FUSION:
```
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
```
* add the bundled GCTA binary `gcta_nr_robust` to path
* download and install [GEMMA software](https://github.com/genetics-statistics/GEMMA/releases) 

* Launch R and install required libraries
```
install.packages(c('optparse','RColorBrewer'))
install.packages('plink2R-master/plink2R/',repos=NULL)
install.packages(c('glmnet','methods'))
```
# Example analysis

## 1st step: Train gene imputation models for each cell type

### Using TIGAR tool to train DPR and Elastic-Net models

- (1) Genotype data: vcf or dosage file
- (2) Training Sample ID file: headerless, single-column file containing sampleIDs to use
- (3) protein abundance file:
-  First 5 columns are _Chromosome number, Gene start position, Gene end position, Target gene ID, Gene name (optional, could be the same as Target gene ID)_

#### Data input
```
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
```

#### Train nonparametric Bayesian Dirichlet Process Regression (DPR) model
```
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
```

#### Train Elastic-Net (PrediXcan) model
```
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
```

### Train FUSION/BestModel model
```
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
```

## 2nd step: Conduct summary-level TWAS using TIGAR for each cell type

#### We integrate cell-type-specific eQTL weights from gene imputation models with GWAS summary statistics to conduct the gene-based association test

- `--asso`: `2` summary-level TWAS using GWAS summary Z-score statistics and reference LD
- `--weight`: Path to SNP weight (eQTL effect size) file
- `--Zscore`: Path to GWAS summary Zscore statistics
- `--LD`: Path to reference LD (SNP genotype covariance matrix) that should be bgzipped and tabixed.
- `--window`: Window size (in base pairs) around gene region from which to include SNPs (default: `1000000` [+- 1MB region around gene region])
- `--test_stat`: burden Z test statistic to calculate: `FUSION`, `SPrediXcan`, or `both` (default both)

```
Zscore=$path_to_gwas_file
LD_file=$path_to_LD_file
TIGAR_dir=$path_to_tigar

weight=$path_to_weight
gene_anno=$path_to_gene_annotation_file
out_dir=$path_to_output

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
```

## 3rd step: Combine p-values from three methods using ACAT-O test
Install the ACAT package in R

```
library(devtools)
devtools::install_github("yaowuliu/ACAT")
```

Launch R and load the package
```
library(ACAT)

ACAT_withNA = function(p_vec){
  p_vec_noNA = p_vec[is.na(p_vec) == F]
  ACAT(p_vec_noNA)
  
}
```

## 4th step: Analyze the CTS-TWAS results
- Generate the Manhattan plot
- Q-Q plot
- Cross-validation `R**2` comparison plot of three gene imputation models

## 5th step: GIFT fine mapping analysis

- (1) eQTL summary statistics for SNPs in a genomic region
- (2) GWAS summary statistics for SNPs in a genomic region
- (3) Reference LD correlation matrix for eQTL
- (4) Reference LD correlation matrix for GWAS

Install the GIFT package in R
```
# Install GIFT package
install.packages('devtools')
devtools::install_github('yuanzhongshang/GIFT')
```

Load the package and run the analysis
```
library(GIFT)
eQTLfilelocation <- $path_to_eQTL
GWASfile <- $path_to_gwas
eQTLLDfile <- $path_to_eQTLLD
GWASLDfile <- $path_to_GWASLD
snplist <- $path_to_snplist
pindex <- $path_to_pindex


convert <- pre_process_summary(eQTLfilelocation, eQTLLDfile, GWASfile, GWASLDfile, snplist, pindex)
gene <- convert$gene
Zscore1 <- convert$Zscore1
Zscore2 <- convert$Zscore2
LDmatrix1 <- convert$LDmatrix1
LDmatrix2 <- convert$LDmatrix2
n1 <- $sample_size_eqtl
n2 <- $sample_size_gwas
R <- $path_to_gene_expressions

result <- GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, n1, n2, gene, pindex, R=R, maxiter=1000, tol=1e-4, pleio=0, ncores=1, in_sample_LD=T)
```

Visualize GIFT results
```
library(ggrepel)
library(ggplot2)
                  
p1 <- ggplot(data) +
  labs(x = "region name", y = expression(paste(-log[10], " (p-value)"))) +
  geom_point(aes(x = BP, y = -log10(P), color = index, shape = index, size = index)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), 
    axis.line.y = element_line(color = "black", linetype = "solid"), 
    axis.line.x = element_line(color = "black", linetype = "solid"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  scale_discrete_manual(values = c("grey", "#377EB8", "#F23557"), aesthetics = 'colour') +
  scale_shape_manual(values = c(19, 15, 18)) +
  scale_size_manual(values = c(1, 1.5, 2)) +
  theme(panel.grid = element_blank()) +
  geom_text_repel(data = subset(data, index == "GIFT" & P < 0.05), 
                  aes(x = BP, y = -log10(P), label = X), 
                  size = 4,
                  fontface="bold",
                  box.padding = 1, 
                  point.padding = 0.8, 
                  segment.color = "black",
                  segment.size = 0.5,
                  nudge_y = 0.3)
```

## Data availability
