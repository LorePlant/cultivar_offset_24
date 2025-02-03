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

plot K2

best = which.min(cross.entropy(pop_stru, K = 2))
qmatrix_K2 = Q(pop_stru, K = 2, run = best)


read.table("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/Pop_list.txt", header = T)
POP_matrixK2<-cbind(POP_info, qmatrix_K2)
POP_matrixK2 <- POP_matrixK2 %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

Create the Bar Plot

K2<-ggplot(POP_matrixK2, aes(x = IND, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) + # Customize colors
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population
theme(
    strip.text.x = element_text(size = 7),  # Reduce POP facet label size
    strip.background = element_blank()  # Remove background for a cleaner look
  )

jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/K2barplot",width = 15, height = 3, units = "cm", res = 800)
K2
dev.off()


plot K3

best = which.min(cross.entropy(pop_stru, K = 3))
qmatrix_K3 = Q(pop_stru, K = 3, run = best)
POP_matrixK3<-cbind(POP_info, qmatrix_K3)
POP_matrixK3 <- POP_matrixK3 %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

Create the Bar Plot

K3<-ggplot(POP_matrixK3, aes(x = IND, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "green")) + # Customize colors
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population
theme(axis.text.x = element_blank(), # Hide individual labels if too many
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) 
  

jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/K3barplot.jpeg",width = 15, height = 3, units = "cm", res = 800)
K3
dev.off()

plot K4

best = which.min(cross.entropy(pop_stru, K = 4))
qmatrix_K4 = Q(pop_stru, K = 4, run = best)
POP_matrixK4<-cbind(POP_info, qmatrix_K4)
POP_matrixK4 <- POP_matrixK4 %>%
  pivot_longer(cols = starts_with("V"), 
               names_to = "Cluster", 
               values_to = "Ancestry")

Create the Bar Plot

K4<-ggplot(POP_matrixK4, aes(x = IND, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Individuals", y = "Ancestry Proportion") +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "green", "darkblue")) + # Customize colors
  facet_grid(~POP, scales = "free_x", space = "free_x") # Separate by population
theme(axis.text.x = element_blank(), # Hide individual labels if too many
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) 
  

jpeg(file = "/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/genotypes.snmf/K4barplot.jpeg",width = 15, height = 3, units = "cm", res = 800)
K4
dev.off()

```











# select the best run for K = 4 clusters

my.colors <- c("tomato", "lightblue",
"olivedrab")
barchart(prova, K = 3, run = best,
border = NA, space = 0,
col = my.colors,
xlab = "Individuals",
ylab = "Ancestry proportions",
main = "Ancestry matrix")

```
