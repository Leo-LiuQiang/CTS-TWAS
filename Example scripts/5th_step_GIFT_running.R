# Install GIFT package
install.packages('devtools')
devtools::install_github('yuanzhongshang/GIFT')

library(GIFT)
#### load the directory containing files of summary statistics from eQTL data only
eQTLfilelocation <- $path_to_eQTL
#### load the directory of summary statistics from GWAS data
GWASfile <- $path_to_gwas
#### load the directory of LD matrix from eQTL data and GWAS data
eQTLLDfile <- $path_to_eQTLLD
GWASLDfile <- $path_to_GWASLD
#### load the SNP list and cis-SNP number for each gene in a region (pindex)
snplist <- $path_to_snplist
pindex <- $path_to_pindex
#### pre-process the file to be a list including gene names vector, z-score matrix and LD matrix of eQTL data and GWAS data
convert <- pre_process_summary(eQTLfilelocation, eQTLLDfile, GWASfile, GWASLDfile, snplist, pindex)
gene <- convert$gene
Zscore1 <- convert$Zscore1
Zscore2 <- convert$Zscore2
LDmatrix1 <- convert$LDmatrix1
LDmatrix2 <- convert$LDmatrix2

### input the sample sizes of eQTL data and GWAS data
n1 <- $sample_size_eqtl
n2 <- $sample_size_gwas

R <- $path_to_gene_expressions

result <- GIFT_summary(Zscore1, Zscore2, LDmatrix1, LDmatrix2, n1, n2, gene, pindex, R=R, maxiter=1000, tol=1e-4, pleio=0, ncores=1, in_sample_LD=T)

# Visualization
#### GWAS result
GWASresult=$path_to_gwas
GWASresult=GWASresult[,c(2,3,11)]
GWASresult$index="GWAS"
colnames(GWASresult)=c("X","BP","P","index")

#### TWAS result
regions <- $path_to_twas
index <- gsub("region", "", region)
regions <- regions[which(regions[,"region"]==index),]
TWASresult=regions[,c(5,1,2,3,10)]
TWASresult$BP=apply(TWASresult[,c(3,4)],1,mean)
TWASresult=TWASresult[,c(1,6,5)]
TWASresult$index="TWAS"
colnames(TWASresult)=c("X","BP","P","index")

#### GIFT result
GIFTresult=result
GIFTresult$BP=TWASresult$BP
GIFTresult=GIFTresult[,c(1,4,3)]
GIFTresult$index="GIFT"
colnames(GIFTresult)=c("X","BP","P","index")

#### visualize the result by Manhattan plot
data=rbind(GWASresult,TWASresult)
data=rbind(data,GIFTresult)
data$BP=data$BP/1000000
data$index=factor(data$index,levels=c("GWAS","TWAS","GIFT"))

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
