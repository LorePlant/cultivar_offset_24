# Landscape Olea europea

This page is created to track progresses on my postdoctoral research in landscape genomics in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.

## vcf file preparation
I started by filtering sites quality

```
vcftools --gzvcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/All_wild_cultivated_olive_2_run.vcf.gz --remove-indels --minDP 10 --max-meanDP 100 --minQ 200  --max-alleles 2 --min-alleles 2 --max-missing 0.90  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf

remouve individual withi missingness >0.85

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf.recode.vcf --missing-indv --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/F_missing_individuals_DP10_100

vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_DP10_100_miss090.vcf.recode.vcf --keep /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/list710WD_cul.txt  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf
```

In total we obtained a dataset of 710 individuals

## Analysis of Population Structure

Let's remouve SNPs that are in linkage using the thinning option in vcf tool. This command will select a SNPs every 1Kbp. Considering the low LD present in Olive (~250bp) we can consider this window appropriate.

```
thinning vcf
vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085.vcf.recode.vcf --thin 1000 --recode --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/WC710_lec24_DP10_100_miss090_ind085_Thinned

```
Analysis of Population Structure using LEA package

```

R
setwd("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24")
geno710 <- read.vcfR("WC710_lec24_DP10_100_miss090_ind085_Thinned.recode.vcf")#import vcf file
GI <- vcfR2genind(geno710)#transfrom file in genind object
geno710<-as.data.frame(gl.genoLAND)
geno710<-geno710%>% select(ends_with(".0"))


write.geno(geno710, "genotypes.geno")

pop_stru = snmf("genotypes.geno", K = 1:15, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/cross_entropy_decay.JPEG")
plot(pop_stru, col = "blue", pch = 19, cex = 1.2)
dev.off()
```
![cross_entropy_decay](https://github.com/user-attachments/assets/a6c19ee9-11bb-4903-bc76-f1d742c207a0)

Print Q matrixes for K runs from K2 to K4 

```

best = which.min(cross.entropy(pop_stru, K = 2))
qmatrix_K2 = Q(pop_stru, K = 2, run = best)

best = which.min(cross.entropy(pop_stru, K = 3))
qmatrix_K3 = Q(pop_stru, K = 3, run = best)


best = which.min(cross.entropy(pop_stru, K = 4))
qmatrix_K4 = Q(pop_stru, K = 4, run = best)

```
Plot bar plot for the K runs using ggplot

```
K2_Qmatrix<-read.table("K2_Qmatrix.txt", header = T)


K2_Qmatrix <- K2_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
  


K2<-ggplot(K2_Qmatrix, aes(x =IND, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("darkblue", "darkorange")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K2

K3_Qmatrix<-read.table("K3_Qmatrix.txt", header = T)


K3_Qmatrix <- K3_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

K3<-ggplot(K3_Qmatrix, aes(x =IND, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("darkblue", "darkorange", "gray")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population

K3

K4_Qmatrix<-read.table("K4_Qmatrix.txt", header = T)


K4_Qmatrix <- K4_Qmatrix%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

K4<-ggplot(K4_Qmatrix, aes(x =IND, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("darkblue", "darkorange", "gray", "darkgreen")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population

K4

ggarrange(K2,K3,K4,nrow=3,ncol=1)

```
![image](https://github.com/user-attachments/assets/735fb6b5-7dfc-41c5-b538-02bdcf8f3e27)

At K4 we can clearly distinguish a specific group present in the Wil dWest and absent in cultivars and Wild East. Using the ancestry coefficient at K=4 we selected individual q>0.7 for the Wild Weast group. 
In total 142 individuals have been selected

We are going to repat the population structure analysis for the 142 individuals with the aim to define population differentiation among pure wild western olive trees.

```
142 WILD WEST
geno142<- read.vcfR("WW142_lec24_DP10_100_miss090_ind085_thinned.vcf.recode.vcf")#import vcf file
G142 <- vcfR2genind(geno142)#transfrom file in genind object
geno142<-as.data.frame(G142)
geno142<-geno142%>% select(ends_with(".0"))
write.geno(geno142, "genotypes.geno")
pop_stru_142WW = snmf("genotypes.geno", K = 1:15, entropy = TRUE, repetitions = 10, project = "new")

# plot cross-entropy criterion for all runs in the snmf project
jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/PopStructure_142WildWest/cross_entropy_decay_142WW.JPEG")
plot(pop_strupop_stru_142WW, col = "blue", pch = 19, cex = 1.2)
dev.off()
```
![image](https://github.com/user-attachments/assets/98646a3a-e471-4f8f-ba03-94b94f8a63e2)


The cross-entropy coefficient reached a minumum at K=5. We are going to use this partition to construnct the ancestry barplot using ggplot and map the spatial interpolation of ancestry coefficients in the species niche using QGIS.

```
best = which.min(cross.entropy(pop_stru_142WW, K = 5))
qmatrix_K5 = Q(pop_stru_142WW, K = 5, run = best)
wild_info<-read.table("142WildWest.txt")
write.table(qmatrix_K5,"/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/PopStructure_142WildWest/QmatrixK2_142_WW.txt")

K5_Q<-read.table("QmatrixK5_142_WW.txt", header = T)
w142_info<-read.table("142W_info.txt", header = T)
K5Q<-cbind(w142_info, K5_Q)

K5Q <- K5Q%>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")
K5W<-ggplot(K5Q, aes(x =IND, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")) + # Customize colors
  theme(
    axis.text.x = element_blank(),  # Hide individual labels if too many
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "lines")) +
  facet_grid(~POP, scales = "free_x", space = "free") # Separate by population

K5W

```
![image](https://github.com/user-attachments/assets/a53b1741-8b0b-45f2-a98c-6131b7ac4921)

![image](https://github.com/user-attachments/assets/cfccd6b5-cf12-4bea-89ae-462a1be31944)


