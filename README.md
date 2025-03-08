
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
test_env <- data_wild%>% select(long, lat, bio2, bio10, bio11, bio15, bio18, bio19, clay, N, pH, sand)
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
RDAgeo_env <- rda(geno155 ~ bio2+bio10+bio11+	bio15	+ bio18 + bio19 + clay+ N+ pH+ sand , Variables)
sqrt(vif.cca(RDAgeo_env))
```

|  bio2    |   bio10   |   bio11  |  bio15  |   bio18  |   bio19  |clay    |   N  |    pH  |     sand |
|---------|----------|---------|----------|-----------|----------|-------|-------|---------|--------|
1.648260 |2.132844| 2.717164| 3.642094| 2.693061 |2.360481|3.067049 |2.399268 |1.567922| 2.776405|


## RDA for Genotype Environment Associations (GEA)

Redundancy analysis can be used to identify GEA based on the Mhallanoise distance of SNPs in the RDA-biplot. Within the RDA model we can effectively correct for population structure  and geography (latitude and longitude) using them as covariates in the RDA model. As population structure correction we used latent factor derived from the LEA package.

As first attempt I decided to run the anlysis seperate for temperature, precipitation and soil variables.

>Temperature
```
Y <- geno155
sel_temp<- data.frame(Env%>% dplyr::select(bio2, bio10, bio11))
write.env(sel_temp, "Temp_variable.env")
X = read.table("Temp_variable.env")

mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3)
str(mod.lfmm2)
mod.lfmm2@U
#Merge latent factor to Variable
latent_temp<-data.frame(rownames(geno155), mod.lfmm2@U)
Temp_Var<-cbind(Variables,latent_temp)

#GEA Temperature
RDA_temp <- rda(geno155 ~ bio2+bio10+bio11 +  Condition(X1 + X2 +X3), Temp_Var)
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
  xlab("RDA 1: 57.8%") + ylab("RDA 2: 21.6%") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_temp
jpeg(file = "/lustre/rocchettil/RDA_temp_biplot.jpeg")
plot(loading_temp)
dev.off()

write.table(qvalue, "Temp_GEA_WWE_3.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

Manhattan_temp <- read.csv(file = "Temp_GEA_WWE_3.csv", header=TRUE) #import the p value result for temperature
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.001158020), genomewideline = -log10(1.924104e-07))
jpeg(file = "/lustre/rocchettil/Manh_RDA_temp.jpeg")
manhattan(Manhattan_temp, col = c("darkred", "gray60"),suggestiveline = -log10(0.001158020), genomewideline = -log10(1.924104e-07))
dev.off()

#P distribution
jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_temp")
hist(Manhattan_temp$P)
dev.off()

hist(qvalue$p.value)
```
![image](https://github.com/user-attachments/assets/19591930-86d0-49a3-8cbf-3c7c78631c31)

GEA for temperature variables resulted in 5877 SNPs  FDR and among them 167 highly assoicciated Bonferroni correction

>Precipitation
  ```
  #latent factor precipitation variable
  
  Y <- geno155
  sel_prec<- data.frame(Env%>% select(bio15, bio18, bio19))
  write.env(sel_prec, "prec_variable.env")
  X = read.table("prec_variable.env")
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3)
  str(mod.lfmm2)
  mod.lfmm2@U
  #Merge latent factor to Variable
  latent_prec<-data.frame(rownames(geno155), mod.lfmm2@U)
  Prec_Var<-cbind(Variables,latent_prec)
  
  
  
  ## GEA Precipitation
  RDA_prec <- rda(geno155 ~ 	bio15	+ bio18 + bio19 +  Condition(X1 + X2 + X3 ), Prec_Var)
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
    xlab("RDA 1: 48.5%") + ylab("RDA 2: 27.0%") +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_prec
  jpeg(file = "/lustre/rocchettil/RDA_prec_biplot.jpeg")
  plot(loading_prec)
  dev.off()
  
  
  write.table(qvalue, "Prec_GEA_WWE_3.csv", append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  #plotting Mhanattan plot using the library qqman
  
  library(qqman)
  Manhattan_prec <- read.csv(file = "Prec_GEA_WWE_3.csv", header=TRUE) #import the p value result for precipitation
  jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
  manhattan(Manhattan_prec, col = c("blue", "gray60"),suggestiveline = -log10(0.0010233235), genomewideline = -log10(1.914289e-07))
  dev.off()
  
  #P distribution
  jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
  hist(Manhattan_prec$P)
  dev.off()

```
![image](https://github.com/user-attachments/assets/e640a7ff-548e-420e-a08b-d535950dff33)

GEA for precipitation variables resulted in 5191 SNPs FDR and among them 177 highly assoicciated Bonferroni correction

>Soil
```
  #latent factor precipitation variable
  
  Y <- geno155
  sel_soil<- data.frame(Env%>% select(clay, N, pH, sand))
  write.env(sel_soil, "soil_variable.env")
  X = read.table("soil_variable.env")
  mod.lfmm2 <- lfmm2(input = Y, env = X, K = 3)
  str(mod.lfmm2)
  mod.lfmm2@U
  #Merge latent factor to Variable
  latent_soil<-data.frame(rownames(geno155), mod.lfmm2@U)
  soil_Var<-cbind(Variables,latent_soil)
  
  
  
  ## GEA Soil
  RDA_soil <- rda(geno155 ~ 	clay+ N+ pH+ sand+ Condition(X1 + X2 + X3 ), soil_Var)
  summary(eigenvals(RDA_soil, model = "constrained"))
  
  rdadapt_soil<- rdadapt(RDA_soil, 2)
  ## P-values threshold after Bonferroni correction
  thres_env <- 0.05/length(rdadapt_soil$p.values)
  ## Identifying the loci that are below the p-value threshold
  top_outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_soil$p.values<thres_env)], p.value = rdadapt_soil$p.values[which(rdadapt_soil$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(geno155)[which(rdadapt_soil$p.values<thres_env)], split = "_"), function(x) x[1])))
  
  qvalue <- data.frame(Loci = colnames(geno155), p.value = rdadapt_soil$p.values, q.value = rdadapt_soil$q.value)
  outliers <- data.frame(Loci = colnames(geno155)[which(rdadapt_soil$q.values<0.05)], p.value = rdadapt_soil$p.values[which(rdadapt_soil$q.values<0.05)])
  
  #plot GEA soil
  
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
    xlab("RDA 1: 39.1%") + ylab("RDA 2: 24.2%") +
    guides(color=guide_legend(title="Locus type")) +
    theme_bw(base_size = 11, base_family = "Times") +
    theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  loading_soil
  jpeg(file = "/lustre/rocchettil/RDA_soil_biplot.jpeg")
  plot(loading_soil)
  dev.off()
  
  
  write.table(qvalue, "Soil_GEA_WWE_3.csv", append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  #plotting Mhanattan plot using the library qqman
  
  library(qqman)
  Manhattan_soil <- read.csv(file = "Soil_GEA_WWE_3.csv", header=TRUE) #import the p value result for precipitation
  jpeg(file = "/lustre/rocchettil/Manh_RDA_prec.jpeg")
  manhattan(Manhattan_soil, col = c("#ab7e4c", "gray60"),suggestiveline = -log10(0.001236238), genomewideline = -log10(1.961360e-07))
  dev.off()
  
  #P distribution
  jpeg(file = "/lustre/rocchettil/Phist_Manh_RDA_prec.jpeg")
  hist(Manhattan_soil$P)
  dev.off()

```
![image](https://github.com/user-attachments/assets/dc0d5d59-148a-4b39-9fa1-ba1076835bff)

GEA for soil variables resulted in 6271 SNPs FDR and among them 194 highly assoicciated Bonferroni correction

## Enriched RDA
To visualize the adaptive differentiation among genotypes, I conducted an additional Redundancy Analysis (RDA) using only the GEA SNPs for the three seperate analysis for temperature, precipitation and soil variables (FDR, q<0.05).
By using the same GEA QTLs, we are going to represent first the RDA biplot with all genotypes including Wild from East and West mediterrenean, than RDA adaptive landscape in west mediterrenean only. This last plot will represent the adaptive landscape in the West mediterrenean that includes potential genetic adaptation derived from easter germplasm
```
 ## A total of 12143 GEA QTL have been identified and used for RDA
  
  geno_Wild_GEA<-geno155[which((rdadapt_temp$q.values<0.05)|(rdadapt_prec$q.values<0.05)|(rdadapt_soil$q.values<0.05))]
  write.table(geno_Wild_GEA, "geno_Wild_GEA_WWE.txt") #save the new GEA genotype data
  geno_Wild_GEA<-read.table("geno_Wild_GEA_WWE.txt")
  
  

  
  RDA_all_enriched<-rda(geno_Wild_GEA ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 + N + pH + clay + sand, Variables)
  summary(eigenvals(RDA_all_enriched, model = "constrained"))
plot(RDA_all_enriched)
sqrt(vif.cca(RDA_all_enriched))

# plot Geographic regions


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

### RDA plot based on GEA considering only West Mediterrenean
#filter west med individuals

list142WW<-read.table("list142WW.txt")
geno_wild_GEA_142WW<- geno_Wild_GEA[rownames(geno_Wild_GEA)%in% list142WW$V1, ]
Variables_142WW<- Variables[Variables$geno%in% list142WW$V1, ]
RDA_142WW_enriched<-rda(geno_wild_GEA_142WW ~ bio2 + bio10 + bio11 + bio15	+ bio18 + bio19 + N + pH + clay + sand, Variables_142WW)
summary(eigenvals(RDA_142WW_enriched, model = "constrained"))

TAB_gen <- data.frame(geno = row.names(scores(RDA_142WW_enriched , display = "sites")), scores(RDA_142WW_enriched, display = "sites"))

Geno <- merge(TAB_gen, Variables[, 1:5] ,by="geno")
TAB_var <- as.data.frame(scores(RDA_142WW_enriched, choices=c(1,2), display="bp"))
loading_geno_142W_enriched_region<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = Geno, aes(x=RDA1, y=RDA2, fill = region), size = 2.5, shape = 21, color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("lightblue","darkgreen", "darkorange")) +
  geom_segment(data = TAB_var, aes(xend=RDA1*5, yend=RDA2*5, x=0, y=0), colour="black", linewidth =0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(data = TAB_var, aes(x=5*RDA1, y=5*RDA2, label = row.names(TAB_var)), size = 3.2, family = "Times") +
  xlab("RDA 1: 41.1%") + ylab("RDA 2: 16.6%") +
  guides(color=guide_legend(title="Latitude gradient")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
loading_geno_142W_enriched_region

ggarrange(loading_geno_all_enriched_region, loading_geno_142W_enriched_region, nrow=1,ncol=2)
```
![image](https://github.com/user-attachments/assets/933ca0b6-e16d-41ac-849a-2baebea0b216)

In the first plot, we can see a clear differentiation among west and east mediterrenean on the first axis driven by soil pH and bio10(summer temperature). On the second axis we see a differentiation among south and north regions of the west mediterrenean, with the norther region in france with higher summer precipitation (bio18)and Nitrogen content opposite to the souther regions (Morocco and Spain) with higher winter temperature and clay.

In the second plot we can see a clear diffentiation between two different landscpae in Morocco potentially highlight the difference betwen continental and costal Morocco. 

## Adaptive landscape projection
By levaraging the enriched RDA model (RDA_142WW_enriched) we can estimate the adaptive value of each pixel within the olive niche.
First lets upload raaster files previously cut using the ENM mask in QGIS
```
library(raster)
library("readxl")


bio2<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio2_current_masked.tif"))
bio10<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio10_current_masked.tif"))
bio11<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio11_current_masked.tif"))
bio15<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio15_current_masked.tif"))
bio18<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio18_current_masked.tif"))
bio19<- raster(paste("D:/raster files/Current_ENM_clipped_biova/bio19_current_masked.tif"))
soilN<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilN.tif"))
soilpH<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilpH.tif"))

soilclay<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilclay.tif"))
soilsand<- raster(paste("D:/raster files/Current_ENM_clipped_biova/resampled_soilsand.tif"))

names(bio2) = 'bio2'
names(bio10) = 'bio10'
names(bio11) = 'bio11'
names(bio15) = 'bio15'
names(bio18) = 'bio18'
names(bio19) = 'bio19'
names(soilN ) = 'N'
names(soilpH) = 'pH'
names(soilclay) = 'clay'
names(soilsand) = 'sand'


#alignment of soil rasters with bioclimatic variables

soilN <- resample(soilN, bio2, method="bilinear")
writeRaster(soilN, "resampled_soilN.tif", format="GTiff", overwrite=TRUE)

soilpH<- resample(soilpH, bio2, method="bilinear")
writeRaster(soilpH, "D:/raster files/Current_ENM_clipped_biova/resampled_soilpH.tif", format="GTiff", overwrite=TRUE)

soilclay <- resample(soilclay, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/Current_ENM_clipped_biova/resampled_soilclay.tif", format="GTiff", overwrite=TRUE)

soilsand <- resample(soilsand, bio2, method="bilinear")
writeRaster(soilclay, "D:/raster files/Current_ENM_clipped_biova/resampled_soilsand.tif", format="GTiff", overwrite=TRUE)

#stack the different raster file
ras_current_var<-stack(c(bio2,bio10, bio11, bio15, bio18, bio19, soilclay,soilN,soilpH, soilsand))
plot(ras_current_var, 
     xlim = c(-10, 12), 
     ylim = c(27, 50))
```
>in this specific code chunk we can find the function _resample_ that I used to allign soil rasters with the bioclim raster

Transform the stacked raster in table and use the previous environmental scaling factor to scale the pixel table.
NB: clay, pH; N and sand were adjusted according to the previous RDA model
```
pixel <- as.data.frame(rasterToPoints(ras_current_var[[row.names(RDA_142WW_enriched$CCA$biplot)]]))
pixel <- data.frame(x=pixel$x, y=pixel$y, bio2=pixel$bio2,bio10=pixel$bio10,bio11=pixel$bio11,bio15=pixel$bio15,bio18=pixel$bio18,bio19=pixel$bio19,clay=pixel$clay/10,N=pixel$N/100,pH=pixel$pH/10,sand=pixel$sand/10)
pixel<-na.omit(pixel)
pixel<- pixel[pixel$x>-10, ]
pixel_env<- pixel%>% dplyr::select(bio2, bio10, bio11, bio15, bio18, bio19,clay, N, pH, sand)

scaled_pixel <- scale(pixel_env, center = env_center, scale = env_scale)
scaled_pixel<-as.data.frame(scaled_pixel)
```
Use the RDA model (_RDA_142WW_enriched_) to predict pixel adaptive value (location within the RDA space)

```

#prediction of pixel in the RDA space
scaled_pixel_LC <- predict(RDA_142WW_enriched, newdata=ENV, type="lc",scaling = "sites")
TAB_pixel_LC<- data.frame(lat = pixel$y, long = pixel$x, scaled_pixel_LC[,1:3])
TAB_var <- as.data.frame(scores(RDA_142WW_enriched, choices=c(1,2), display="bp"))
```
Plot the biplot and the geographic projection

```
# Extract RDA values
a1 <- TAB_pixel_LC$RDA1
a2 <- TAB_pixel_LC$RDA2

# Compute the distance from the origin
distance <- sqrt(a1^2 + a2^2)

# Assign colors based on quadrants and the 5th sector (circle radius < 0.5)
TAB_pixel_LC$color <- ifelse(distance < 0.5, "#717171",  # 5th sector - Purple
                             ifelse(a1 > 0 & a2 > 0, "#E41A1C",  # Quadrant 1 - Red
                                    ifelse(a1 < 0 & a2 > 0, "#377EB8",  # Quadrant 2 - Blue
                                           ifelse(a1 < 0 & a2 < 0, "#4DAF4A",  # Quadrant 3 - Green
                                                  "#FF7F00"))))  # Quadrant 4 - Orange

# Update ggplot with quadrant-based colors and 5th sector
pp <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = gray(0.80), size = 0.6) +
  geom_point(data = TAB_pixel_LC, aes(x = RDA1, y = RDA2, color = color), size = 2) +  # Use sector colors
  geom_segment(data = TAB_var, aes(xend = RDA1*3, yend = RDA2*3, x = 0, y = 0), 
               colour = "black", size = 0.15, linetype = 1, 
               arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
  geom_label_repel(data = TAB_var, aes(x = RDA1*3, y = RDA2*3, label = row.names(TAB_var)), 
                   size = 4, family = "Times") +
  xlab("RDA 1: 41.1%") + 
  ylab("RDA 2: 16.6%") +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text = element_text(size = rel(0.8)), 
        strip.text = element_text(size = 9)) +
  scale_color_identity()  # Use predefined colors directly

pp
## plot in geographic map

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Load geographic boundaries of France, Spain, Morocco, Portugal, and Algeria
countries <- ne_countries(scale = "medium", country = c("France", "Spain", "Morocco", "Portugal", "Algeria"), returnclass = "sf")

# Remove French Guiana and Atlantic French territories
countries <- countries[!(countries$geounit %in% c("French Guiana", "Guadeloupe", "Martinique", "Saint Pierre and Miquelon", 
                                                  "Reunion", "Mayotte", "New Caledonia", "French Polynesia", 
                                                  "Wallis and Futuna", "Saint Barthelemy", "Saint Martin")), ]

# Convert TAB_pixel_LC to an sf object
TAB_pixel_LC_sf <- st_as_sf(TAB_pixel_LC, coords = c("long", "lat"), crs = 4326)  # Assumes 'longitude' and 'latitude' columns exist
# Create the map
map <- ggplot(data = countries) +
  geom_sf(fill = "#EBEBEB", color = "black") +  # Countries' borders
  geom_sf(data = TAB_pixel_LC_sf, aes(color = color), size = 0.05, show.legend = FALSE) +  # Points with custom colors
  scale_color_identity() +  # Use exact colors from the 'color' column
  coord_sf(xlim = c(-15, 15), ylim = c(28, 52), expand = FALSE) +  # Set geographic limits
  theme_minimal() +
  labs(title = "Adaptive Landscape") +
  theme(panel.background = element_blank())

map
```
![image](https://github.com/user-attachments/assets/3565dc05-4468-44fc-ba7b-b628c9565158)
