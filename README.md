# Landscape Olea europea

This page is created to track progresses on my postdoctoral research in landscape genomics in a wester Mediterrenean Olive population.
The population is composed by 359 individuals along a 15° latitude gradient from 30 to 45.

## vcf file preparation
I started by filtering sites quality

```
vcftools --gzvcf /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/All_wild_cultivated_olive_2_run.vcf.gz --remove-indels --minDP 8 --max-meanDP 400 --minQ 200  --max-alleles 2 --min-alleles 2 --max-missing 0.90  --recode --recode-INFO-all --out /storage/replicated/cirad/projects/CLIMOLIVEMED/results/GenomicOffsets/Lorenzo/Leccino_new_genome24/site_filtered_wild_cult_lec24_miss_090.vcf
```

Filter individual with high missngness 
