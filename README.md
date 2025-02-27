
This page is created to track progresses on my postdoctoral research in landscape genomics in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.

## vcf file preparation
I started by filtering sites quality

```
vcftools --gzvcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/All_wild_cultivated_olive_2_run.vcf.gz --remove-indels --minDP 10 --max-meanDP 100 --minQ 200  --max-alleles 2 --min-alleles 2 --max-missing 0.90  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf

remouve individual withi missingness >0.85

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf.recode.vcf --missing-indv --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/F_missing_individuals_DP10_100

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf.recode.vcf --keep /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/list710WD_cul.txt  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf

further remouve OES_M29_07_S9_L004  OES_M29_02_S56_L004

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf --keep /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/list710WD_cul.txt  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC708_lec24_DP10_100_miss090_ind085.vcf

Found monomorphic SNPs. Filter for mac minor allele count 1

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC708_lec24_DP10_100_miss090_ind085.vcf --mac 1  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC708_lec24_DP10_100_miss090_ind085_mac1.vcf

```

In total we obtained a dataset of 708 individuals

## Analysis of Population Structure

Let's remouve SNPs that are in linkage using the thinning option in vcf tool. This command will select a SNPs every 1Kbp. Considering the low LD present in Olive (~250bp) we can consider this window appropriate.

```
thinning vcf
vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf.recode.vcf --thin 1000 --recode --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085_Thinned

```
Analysis of Population Structure using LEA package

```
geno708 <- read.vcfR("WC708_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf")#import vcf file
GI <- vcfR2genind(geno708)#transfrom file in genind object
geno708<-as.data.frame(GI)
geno708<-geno708%>% select(ends_with(".0"))
list708<-data.frame(row.names(geno708))
write.table(list708, "list708.txt")#save individual order

write.geno(geno708, "Pop_stru_708.geno")

pop_stru = snmf("Pop_stru_708.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")
```
![image](https://github.com/user-attachments/assets/2c7f1026-c05f-4253-8d50-b910331a0c9b)


Print Q matrixes for K runs from K2 to K4 

```

best = which.min(cross.entropy(pop_stru, K = 2))
qmatrix_K2 = Q(pop_stru, K = 2, run = best)

best = which.min(cross.entropy(pop_stru, K = 3))
qmatrix_K3 = Q(pop_stru, K = 3, run = best)


best = which.min(cross.entropy(pop_stru, K = 4))
qmatrix_K4 = Q(pop_stru, K = 4, run = best)


pop_info_708<-read.table("708_pop_info.txt", header = T)

```
Plot bar plot for the K runs using ggplot

```
K2_Qmatrix<-cbind(pop_info_708, qmatrix_K2)


K2_Qmatrix <- K2_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
  


K2<-ggplot(K2_Qmatrix, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("deepskyblue4", "darkorange")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K2

K3_Qmatrix<-cbind(pop_info_708, qmatrix_K3)


K3_Qmatrix <- K3_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

K3<-ggplot(K3_Qmatrix, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("darkorange","deepskyblue4", "darkgray")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population

K3

K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)


K4_Qmatrix <- K4_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

K4<-ggplot(K4_Qmatrix, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("darkgreen", "darkorange", "gray", "deepskyblue4")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population

K4

ggarrange(K2,K3,K4,nrow=3,ncol=1)

```
![image](https://github.com/user-attachments/assets/22f6cad6-556c-41ea-a115-40603761bd78)


At K4 we can clearly distinguish a specific group present in the Wil West and absent in cultivars and Wild East. Using the ancestry coefficient at K=4 we selected individual q>0.7 for the Wild Weast group. 
In total 142 individuals have been selected

```
K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)
pure_wild_west <- subset(K4_Qmatrix, V1  > 0.7)
write.table(pure_wild_west, "pure_wildW_070.txt")
```

We are going to repat the population structure analysis for the 142 individuals with the aim to define population differentiation among pure wild western olive trees.

```
# pure wild based on K=4


Pop_stru_708 <- load.snmfProject("Pop_stru_708.snmfProject")

best = which.min(cross.entropy(Pop_stru_708, K = 4))
qmatrix_K4 = Q(Pop_stru_708, K = 4, run = best)
K4_Qmatrix<-cbind(pop_info_708, qmatrix_K4)
pure_wild_west <- subset(K4_Qmatrix, V1  > 0.7)
pure_wildW <- pure_wild_west %>% select(id)

geno142_WW <-  geno708[rownames(geno708)%in% pure_wildW$id, ]
list142_wildW<- data.frame(rownames(geno142_WW))
write.table(list142_wildW, "list_142WW.txt")
setwd("C:/Users/rocchetti/Desktop/Leccino24/PopulationStructure/Pop_structure_142_wild")

write.geno(geno142_WW, "Pop_stru_142_WW.geno")
pop_stru_142WW = snmf("Pop_stru_142_WW.geno", K = 1:10, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/cross_entropy_decay.JPEG")
plot(pop_stru_142WW, col = "blue", pch = 19, cex = 1.2)
dev.off()

```
![image](https://github.com/user-attachments/assets/44d90bd8-540d-4a1d-959b-97dc463c5223)




The cross-entropy coefficient reached a minumum at K=3. We are going to use this partition to construnct the ancestry barplot using ggplot and map the spatial interpolation of ancestry coefficients in the species niche using QGIS.

```
best = which.min(cross.entropy(pop_stru_142WW, K = 2))
K2Q = Q(pop_stru_142WW, K = 2, run = best)
write.table(K2Q, "Qmatrix_K2_142WW.txt")

best = which.min(cross.entropy(pop_stru_142WW, K = 3))
K3Q = Q(pop_stru_142WW, K = 3, run = best)
write.table(K3Q, "Qmatrix_K3_142WW.txt")

pop_info_142WW<- pop_info_708[pop_info_708$id %in% pure_wildW$id, ]

K2Q<-cbind(pop_info_142WW, K2Q)

K2Q <- K2Q%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
K2W<-ggplot(K2Q, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K2W


K3Q<-cbind(pop_info_142WW, K3Q)

K3Q <- K3Q%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
K3W<-ggplot(K3Q, aes(x =id, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "darkgrey")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K3W

```
![image](https://github.com/user-attachments/assets/16fb6a37-bc28-4731-9270-a2f64100a53c)

![image](https://github.com/user-attachments/assets/255ddb7e-5b89-4247-ad17-0003e53ecb8a)


## hybrid index

142 truly wild genotypes were selected with ancestry q>0.70. To distinguish between historical vs recent introgression, we analyzed the hybrid index using ancestry-informative SNPs. We run the intercalls heterozygosity analysis classyfing as parental group the two wild gene pool Wild East and Wild West.
The significant admixture present in our collection between wild and cultivated material can be derived from recent crossing forming F1 hybrids or from past generations of crossing where natural selection had the possibilty to act. The GEA identification presume that the associated QTL derived from processes of local adaptation where environmental selection had the generational time to act. In our case the cultivated genome from cultivars vegetatevly progated mainly evolved in the eastern part of the mediterrenan basin, so it is paramount to identify the presence of recent F1 hybrids where the cultivated genome have been recently introduced and where selections did not have the generational time to act.

To investigate the presence of F1 hybrids I identified a recent devoped Rpackage that allow to identified ancestry-informative markers and estimate their hybrids index with the relative presence of F1, BC1, BC2 or past introgression. https://omys-omics.github.io/triangulaR/index.html

In this analysis I used the vcf file that was not filtered for MAF 5%.

```
library(triangulaR)
# make a pop map
popmap<-read.table("pop_map_wild_adm.txt", header = T)

  genoLAND.VCF <- read.vcfR("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085_Thinned.recode.vcf")#import vcf file

library(triangulaR)

setwd("C:/Users/rocchetti/Desktop/Leccino24/hybridization_Wild_cult")
popmap<-read.table("pop_map_wild_adm.txt", header = T)
geno708 <- read.vcfR("WC708_lec24_DP10_100_miss090_ind085_mac1_Thinned.recode.vcf")#import vcf file
vcfR.diff <- alleleFreqDiff(vcfR = geno708, pm = popmap, p1 = "WW", p2 = "WE", difference = 0.8)
# 1405 sites passed allele freauency differences of 0.8

hi.het <- hybridIndex(vcfR = vcfR.diff, pm = popmap, p1 = "WW", p2 = "WE")
cols <- c("darkgrey", "purple", "darkorange", "darkgreen")
triangle.plot(hi.het, colors = cols)
```
![image](https://github.com/user-attachments/assets/bcb97c23-7c9f-463f-8bae-67de89a23833)

The results highlight the large presence of rencet hybrids like F1 and BC1, making the admixed population not suited for landscape genomics and GEA discovery. The potential adaptation of these individuals can be given from phenotypic plasticity and/or hybrid vigor.
The large part of cultivated material is closer to the WildEast group. Among them only a few hybrids, mainly from Italt Spain and France (Picholine) are considered F1 with the WildWest group. These result confirm the unexplored diversity of WildWest for the Cultivated germplasm.
Similar results are confirmed as well from PCA analysis

```
#PCA
library(FactoMineR)
library(factoextra)

res.pca708<-PCA(geno708, scale.unit = TRUE, ncp = 5, graph = TRUE)
ind708 <- get_pca_ind(res.pca708)
pca_data708 <- as.data.frame(ind708$coord)
pca_data708<-cbind(popmap, pca_data708)
qq<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = pca_data708, aes(x=Dim.1, y=Dim.2, color = pop), size = 2.5) +
  scale_color_manual(values = c("darkgrey", "purple", "darkorange", "darkgreen")) +
  xlab("PC1: 10%") + ylab("PC2: 5.9%") +
  guides(color=guide_legend(title="Group")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
qq
```
![image](https://github.com/user-attachments/assets/581a716a-2c28-46c4-a84c-35a2cab53d3f)

## RDA

filtered for MAF

R code
```
#enter vcf file
geno155 <- read.vcfR("D:/vcf_file_GEA_leccino/WC156_lec24_DP10_100_miss090_ind085_mac1_MAF005.vcf.recode.vcf")#import vcf file
GI <- vcfR2genind(geno155)#transfrom file in genind object
geno155<-as.data.frame(GI)
geno155<-geno155%>% select(ends_with(".0"))
#imputation
for (i in 1:ncol(geno155))
{
  geno155[which(is.na(geno155[,i])),i] <- median(geno155[-which(is.na(geno155[,i])),i], na.rm=TRUE)
}
geno155_data<- write.table(geno155, "geno_155.txt")
# Wild Environment datafile

#standardize bioclim variable
data_wild<- read.csv("Env_155_WWE.csv", header = TRUE)
test_env <- data_wild%>% select(long, lat, bio2, bio10, bio11, bio15, bio18, bio19)
Env <- scale(test_env, center=TRUE, scale=TRUE)
# Extract the centering values
env_center <- attr(Env, "scaled:center") #mean of each variable
# Extract the scaling values
env_scale <- attr(Env, "scaled:scale") #standard deviation of each variable
#transform into dataset
Env <- as.data.frame(Env)


#combining geographic, Popstructure, environmental (scaled) variables
Variables <- data.frame(data_wild$id, data_wild$group,data_wild$region, Env)
names(Variables)[1]<- paste("geno")
names(Variables)[2]<- paste("group")
names(Variables)[3]<- paste("region")
```
Check Variance Inflation Factor of selected environmental variable including bioclim and soil variable;

```
RDAgeo_env <- rda(genotype ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19, Variables)

sqrt(vif.cca(RDAgeo_env))
```

|  bio2    |   bio10   |   bio11  |  bio15  |   bio18  |   bio19  | N    |   pH  |    SOC  |     Wc |
|---------|----------|---------|----------|-----------|----------|-------|-------|---------|--------|
2.264430 |1.849025| 2.164724| 5.224877| 3.740733 |2.981721|3.402853 |1.710886 |4.158232| 2.253560|


## RDA for Genotype Environment Associations (GEA)

Redundancy analysis can be used to identify GEA based on the Mhallanoise distance of SNPs in the RDA-biplot. Within the RDA model we can effectively correct for population structure  and geography (latitude and longitude) using them as covariates in the RDA model. As population structure correction we used latent factor derived from the LEA package.

As first attempt I decided to run the anlysis seperate for temperature, precipitation and soil variables.

>Temperature
```
## Use latent factor for covariable correction
# latent factor temperature variable
Y <- geno155
sel_temp<- data.frame(Env%>% dplyr::select(bio2, bio10, bio11))
write.env(sel_temp, "Temp_variable.env")
X = read.table("Temp_variable.env")

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 4)
str(mod.lfmm2)
mod.lfmm2@U
#Merge latent factor to Variable
latent_temp<-data.frame(rownames(geno155), mod.lfmm2@U)
Temp_Var<-cbind(Variables,latent_temp)

#GEA Temperature
RDA_temp <- rda(geno155 ~ bio2+bio10+bio11 +  Condition(X1 + X2 +X3 +X4), Temp_Var)
summary(eigenvals(RDA_temp, model = "constrained"))
library(robust)
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

rdadapt_temp<- rdadapt(RDA_temp, 2)
## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_temp$p.values)
## Identifying the loci that are below the p-value threshold
top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_temp$p.values<thres_env)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_temp$p.values<thres_env)], split = "_"), function(x) x[1])))
qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_temp$p.values, q.value = rdadapt_temp$q.value)
outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_temp$q.values<0.05)], p.value = rdadapt_temp$p.values[which(rdadapt_temp$q.values<0.05)])


#plot GEA temp

locus_scores <- scores(RDA_temp, choices=c(1:2), display="species", scaling="sites")
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Not associated"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
TAB_var <- as.data.frame(scores(RDA_temp, choices=c(1,2), display="bp"))
loading_temp<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1, y=RDA2, colour = type), size = 2.5) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=1.1*RDA1, yend=1.1*RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 3.8, family = "Times") +
  xlab("RDA 1: 58.2%") + ylab("RDA 2: 21.6%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_temp
jpeg(file = "/lustre/rocchettil/RDA_temp_biplot.jpeg")
plot(loading_temp)
dev.off()

write.table(qvalue, "Temp_GEA_WWE.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

Manhattan_temp <- read.csv(file = "Temp_GEA_WWE.csv", header=TRUE) #import the p value result for temperature
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.001224640), genomewideline = -log10(1.924104e-07))
jpeg(file = "/lustre/rocchettil/Manh_RDA_temp.jpeg")
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.001224640), genomewideline = -log10(1.924104e-07))
dev.off()

#P distribution
jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_temp")
hist(Manhattan_temp$P)
dev.off()

hist(qvalue$p.value)
```
![image](https://github.com/user-attachments/assets/a32731fa-8f4a-4cc6-8596-318219723a26)
GEA for temperature variables resulted in 6213 SNPs  FDR and among them 172 highly assoicciated Bonferroni correction

>Precipitation
  ```
  #latent factor precipitation variable
  
  Y <- geno155
  sel_prec<- data.frame(Env%>% select(bio15, bio18, bio19))
  write.env(sel_prec, "prec_variable.env")
  X = read.table("prec_variable.env")
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = 4)
  str(mod.lfmm2)
  mod.lfmm2@U
  #Merge latent factor to Variable
  latent_prec<-data.frame(rownames(geno155), mod.lfmm2@U)
  Prec_Var<-cbind(Variables,latent_prec)
  
  ## GEA Precipitation
  RDA_prec <- rda(geno155 ~ 	bio15	+ bio18 + bio19 +  Condition(X1 + X2 + X3 + X4 ), Prec_Var)
  summary(eigenvals(RDA_prec, model = "constrained"))
  
  rdadapt_prec<- rdadapt(RDA_prec, 2)
  ## P-values threshold after Bonferroni correction
  thres_env <- 0.05/length(rdadapt_prec$p.values)
  ## Identifying the loci that are below the p-value threshold
  top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_prec$p.values<thres_env)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_prec$p.values<thres_env)], split = "_"), function(x) x[1])))
  
  qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_prec$p.values, q.value = rdadapt_prec$q.value)
  outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_prec$q.values<0.05)], p.value = rdadapt_prec$p.values[which(rdadapt_prec$q.values<0.05)])
  
  #plot GEA precipitation
  locus_scores <- scores(RDA_prec, choices=c(1:2), display="species", scaling="sites")
  TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
  TAB_loci$type <- "Not associated"
  TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
  TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
  TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
  TAB_var <- as.data.frame(scores(RDA_prec, choices=c(1,2), display="bp"))
  loading_prec<-ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_point(data = TAB_loci, aes(x=RDA1, y=RDA2, colour = type), size = 2.5) +
    scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
    geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.8, family = "Times") +
    xlab("RDA 1: 47.5%") + ylab("RDA 2: 27.6%") +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_prec
  jpeg(file = "/lustre/rocchettil/RDA_prec_biplot.jpeg")
  plot(loading_prec)
  dev.off()

  write.table(qvalue, "Prec_GEA_WWE.csv", append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  #plotting Mhanattan plot using the library qqman
  
  library(qqman)
  Manhattan_prec <- read.csv(file = "Prec_GEA_WWE.csv", header=TRUE) #import the p value result for precipitation
  jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
  manhattan(Manhattan_prec, col = c("blue", "gray60"),suggestiveline = -log10(0.0010843300), genomewideline = -log10(1.914289e-07))
  dev.off()
  
  #P distribution
  jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
  hist(Manhattan_prec$P)
  dev.off()
```
![image](https://github.com/user-attachments/assets/d11b4e72-5c10-489f-b12b-2c7dc7ef6b1c)
GEA for precipitation variables resulted in 5501 SNPs FDR and among them 182 highly assoicciated Bonferroni correction

>Soil
```
  #latent factor precipitation variable
  
  Y <- geno155
  sel_soil<- data.frame(Env%>% select(N,pH,SOC,Wc))
  write.env(sel_soil, "soil_variable.env")
  X = read.table("soil_variable.env")
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = 4)
  str(mod.lfmm2)
  mod.lfmm2@U
  #Merge latent factor to Variable
  latent_soil<-data.frame(rownames(geno155), mod.lfmm2@U)
  soil_Var<-cbind(Variables,latent_soil)
  
  
  
  ## GEA Precipitation
  RDA_soil <- rda(geno155 ~ 	N +pH + SOC + Wc +  Condition(X1 + X2 + X3 + X4 ), soil_Var)
  summary(eigenvals(RDA_soil, model = "constrained"))
  
  rdadapt_soil<- rdadapt(RDA_soil, 2)
  ## P-values threshold after Bonferroni correction
  thres_env <- 0.05/length(rdadapt_soil$p.values)
  ## Identifying the loci that are below the p-value threshold
  top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_soil$p.values<thres_env)], p.value = rdadapt_soil$p.values[which(rdadapt_soil$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_soil$p.values<thres_env)], split = "_"), function(x) x[1])))
  
  qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_soil$p.values, q.value = rdadapt_soil$q.value)
  outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_soil$q.values<0.05)], p.value = rdadapt_soil$p.values[which(rdadapt_soil$q.values<0.05)])
  
  #plot GEA precipitation
  
  locus_scores <- scores(RDA_soil, choices=c(1:2), display="species", scaling="sites")
  TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
  TAB_loci$type <- "Not associated"
  TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "FDR"
  TAB_loci$type[TAB_loci$names%in%top_outliers$Loci] <- "Bonferroni"
  TAB_loci$type <- factor(TAB_loci$type, levels = c("Not associated", "FDR", "Bonferroni"))
  TAB_var <- as.data.frame(scores(RDA_soil, choices=c(1,2), display="bp"))
  loading_soil<-ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
    geom_point(data = TAB_loci, aes(x=RDA1, y=RDA2, colour = type), size = 2.5) +
    scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
    geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
    geom_label_repel(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 3.8, family = "Times") +
    xlab("RDA 1: 40.0%") + ylab("RDA 2: 26.8%") +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_soil
  jpeg(file = "/lustre/rocchettil/RDA_soil_biplot.jpeg")
  plot(loading_soil)
  dev.off()
  
  
  write.table(qvalue, "Soil_GEA_WWE.csv", append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  #plotting Mhanattan plot using the library qqman
  
  library(qqman)
  Manhattan_soil <- read.csv(file = "Soil_GEA_WWE.csv", header=TRUE) #import the p value result for precipitation
  jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
  manhattan(Manhattan_soil, col = c("#ab7e4c", "gray60"),suggestiveline = -log10(0.001560509), genomewideline = -log10(1.961360e-07))
  dev.off()
  
  #P distribution
  jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
  hist(Manhattan_soil$P)
  dev.off()
```
![image](https://github.com/user-attachments/assets/8c6d0dc2-0e05-4991-8b2c-fa4cfd69948b)
GEA for soil variables resulted in 7916 SNPs FDR and among them 287 highly assoicciated Bonferroni correction

## Enriched RDA
To visualize the adaptive differentiation among genotypes, I conducted an additional Redundancy Analysis (RDA) using only the GEA SNPs for the three seperate analysis for temperature, precipitation and soil variables (FDR, q<0.05).
```
  geno_Wild_GEA<-geno155[which((rdadapt_temp$q.values<0.05)|(rdadapt_prec$q.values<0.05)|(rdadapt_soil$q.values<0.05))]
  write.table(geno_Wild_GEA, "geno_Wild_GEA_WWE.txt") #save the new GEA genotype data
```
 A total of 13856 GEA QTLs where identified combining Temeprature, precipitation and soil variables (FDR q<0.05). 
This set of GEA QTL will be used as genomic information for running enriched RDA
```
 RDA_all_enriched<-rda(geno_Wild_GEA ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19, Variables)
  summary(eigenvals(RDA_all_enriched, model = "constrained"))
plot(RDA_all_enriched)
```
We re-arrange the plot to highlight specific geographic region
```
TAB_gen <- data.frame(geno = row.names(scores(RDA_all_enriched , display = "sites")), scores(RDA_all_enriched, display = "sites"))

Geno <- merge(TAB_gen, Variables[, 1:5] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_all_enriched, choices=c(1,2), display="bp"))
loading_geno_all_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = region), size = 2.5, shape = 21, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("lightblue","darkgreen", "darkorange", "darkblue")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
  xlab("RDA 1: 40.2%") + ylab("RDA 2: 32.2%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_all_enriched_region
jpeg(file = "/lustre/rocchettil/RDA_all_geno_biplot_region.jpeg")
plot(loading_geno_all_enriched_region)
dev.off()
```
![image](https://github.com/user-attachments/assets/36cec944-f628-4ab7-9f2e-614458230edc)

we can see a clear differentiation among west and east mediterrenean on the first axis driven by soil pH and bio10(summer temperature). On the second axis we see a differentiation among south and noth region, with the norther region in france with higher summer precipitation (bio18), Nitrogen content, organinc matter and soil water capacity.




