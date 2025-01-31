# Landscape Olea europea

This page is created to track progresses on my postdoctoral research in landscape genomics in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.

## vcf file preparation
I started by filtering sites quality

```
vcftools --gzvcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/All_wild_cultivated_olive_2_run.vcf.gz --remove-indels --minDP 8 --max-meanDP 400 --minQ 200  --max-alleles 2 --min-alleles 2 --max-missing 0.90  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_miss_090.vcf
```

Filter individual with high missngness 


Tentative LEA package



```
thinning vcf
vcftools --vcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_miss_095.vcf.recode.vcf --thin 1000 --recode --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/thinned_prova
R
setwd("/storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24")
geno769 <- read.vcfR("thinned_prova.recode.vcf")#import vcf file
gl.genoLAND <- vcfR2genind(geno769)#transfrom file in genind object
geno769<-as.data.frame(gl.genoLAND)
geno769<-geno769%>% select(ends_with(".0"))


write.geno(geno769, "genotypes.geno")

prova = snmf("genotypes.geno",
K = 1:3,
entropy = TRUE,
repetitions = 10,
project = "new")

qmatrix = Q(prova, K = 3, run = best)



# select the best run for K = 4 clusters
best = which.min(cross.entropy(prova, K = 3))
my.colors <- c("tomato", "lightblue",
"olivedrab")
barchart(prova, K = 3, run = best,
border = NA, space = 0,
col = my.colors,
xlab = "Individuals",
ylab = "Ancestry proportions",
main = "Ancestry matrix")

```
