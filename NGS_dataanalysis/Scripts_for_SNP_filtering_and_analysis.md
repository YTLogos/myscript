
SNP filtering for all datasets, 23 May 2017 

### remove indels
```
vcftools --gzvcf MiHiNew_NewOnly2.recode.vcf.gz --remove-indels --recode --recode-INFO-all --out MiHiNew_NewOnly2_noindels
gzip MiHiNew_NewOnly2_noindels.recode.vcf
```
### select only bialelic sites
```
awk '$5 !~ /([[:alpha:]])+,[[:alpha:]]/{print}' MiHiNew_NewOnly2_noindels.recode.vcf > MiHiNew_NewOnly2_noindels_biallelic.vcf
gzip MiHiNew_NewOnly2_noindels_biallelic.vcf
```
### calculate missingness
```
vcftools --gzvcf MiHiNew_NewOnly2_noindels_biallelic.vcf.gz \
--missing-indv \
--out MiHiNew_NewOnly2_noindels_biallelic.vcf
```
### Create list of sampes that meet missingness criteria (from that calculated above)
```
awk ' $5 < 0.9 ' MiHiNew_NewOnly2_noindels_biallelic.vcf.imiss > MiHiNew_NewOnly2_noindels_biallelic_less90miss
awk '{ print $1 }' MiHiNew_NewOnly2_noindels_biallelic_less90miss > List_MiHiNew_NewOnly2_noindels_biallelic_less90miss
```
### Filter them out
```
vcftools --gzvcf MiHiNew_NewOnly2_noindels_biallelic.vcf.gz \
--keep List_MiHiNew_NewOnly2_noindels_biallelic_less90miss \
--recode \
--recode-INFO-all \
--out MiHiNew_NewOnly2_noindels_biallelic_less90miss
gzip MiHiNew_NewOnly2_noindels_biallelic_less90miss.recode.vcf
```

### filter by min/max allele freq cut-off (95 and 5)
```
vcftools --gzvcf MiHiNew_NewOnly2_noindels_biallelic_less90miss.recode.vcf.gz \
--maf 0.05 \
--max-maf 0.95 \
--recode \
--recode-INFO-all \
--out MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95
```

### Apply a "Hard" filter to the dataset
```
bgzip MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95.recode.vcf
tabix MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95.recode.vcf.gz

java -Xmx12G -jar GATK-3.3/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R pgl_assembly_v1.fasta \
-V MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95.recode.vcf.gz \
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 60.0" \
--filterName "FS" \
--filterExpression "MQ < 40.0" \
--filterName "MQ" \
--filterExpression "MQRankSum < -12.5" \
--filterName "MQRankSum" \
--filterExpression "ReadPosRankSum < -8.0" \
--filterName "ReadPosRankSum" \
-o MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfilterannotated.vcf
```

### apply the hard fitler - remove SNPs that fail
```
awk '/^#/||$7=="PASS"' MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfilterannotated.vcf > MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered.vcf
```

### Remove variants with a GQ < 20
```
vcftools --vcf MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered.vcf \
--minGQ 20 \
--recode \
--recode-INFO-all \
--out MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20
```



### Filtering for RAD1 dataset and subsequent analyses
### (Species-level filtering for comparing pop1 and pop2 vs pop7 and pop8)

### filter to just species before selecting SNPs based on cov criteria
```
vcftools --vcf MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20.recode.vcf \
--keep PgPc_groupedasspecies \
--recode \
--recode-INFO-all \
--out MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly
```
### first bgzip then index file
```
set path = ( $path /afs/crc.nd.edu/group/hellmann/hlm_3/Seans_Genomics_data_PgPc/tabix-0.2.6/tabix-0.2.6 )
bgzip MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly.recode.vcf
tabix MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly.recode.vcf.gz
```

### find sites with 75% having 6X cov
```
module load java/1.7
java -Xmx12G -jar GATK-3.3/GenomeAnalysisTK.jar \
-T CoveredByNSamplesSites \
-R pgl_assembly_v1.fasta \
--variant MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly.recode.vcf.gz \
-minCov 6 \
--percentageOfSamples .75 \
-out MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind
-nt 6
```
### need to make tab delimited
```
awk -F":" '$1=$1' OFS="\t" MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind > SNPs_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind.txt

### now filter for just these sites
vcftools --gzvcf MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly.recode.vcf.gz \
--positions SNPs_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind.txt \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind
```

### remove sites that are not polymorphic (for some reason PGDspider is not doing this); can use the 3X because it includes the 6X
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind.recode.vcf \
--exclude-positions sites_monomorphic_3X \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly
```
### CONVERT VCF TO BAYESCAN FILE FORMAT
```
module load java/1.7
java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar \
-inputfile Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly.recode.vcf \
-inputformat VCF \
-outputfile Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly_2pops.bayescanfrmt \
-outputformat GESTE_BAYE_SCAN \
-spid VCFtoBAYESCAN.spid
```
### CONVERT VCF TO ARLEQUIN FILE FORMAT
```
module load java/1.7
java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar \
-inputfile Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly.recode.vcf \
-inputformat VCF \
-outputfile Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly_2pops.arp \
-outputformat ARLEQUIN \
-spid VCFtoARLEQUIN.spid
```

```
### need to add this to the end of the Arlequin file
#
#
#
#        [[Structure]]
#
#                StructureName="New Edited Structure"
#                NbGroups=2
#
#                Group={
#                        "1"
#                        "2"
#                }
#
#                Group={
#                        "7"
#                        "8"
#                }



```

### DETERMINE CHROMOSOME ASSIGNMENT
```
#make dictionary file for protein database
makeblastdb -in silkpep.fa -dbtype prot
#blast Pg to Bmori using all protein sequences of both (whole genome)
./blastp -db silkpep.fa -out All_Bmori_prot_to_All_Pg_prot -query glaucus_official_proteins.fa -evalue 1e-10 -outfmt 6
#now select just the first hit (most sig.) and print it (MUST BE IN BASH TO DO THIS)
awk '!a[$1]++' All_Bmori_prot_to_All_Pg_prot > All_Bmori_prot_to_All_Pg_prot_1stHitOnly
#use the protein fasta of Bmori to create a lookup table to determine nscaf for best mapped protein seqs
grep '>' silkpep.fa > silkpep_used_for_nscaf_info.txt
```




### calculate Fst
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly.recode.vcf \
--weir-fst-pop Pg_species_list \
--weir-fst-pop Pc_species_list \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly_PgvsPc
```

### calculate pi for (pure) P. glaucus populations 
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly.recode.vcf \
--keep Pg_species_list \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly_Pg
```
### calculate pi for (pure) P. canadensis populations 
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly.recode.vcf \
--keep Pc_species_list \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly_Pc
```

### calculate het for all individuals
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly.recode.vcf \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_PgPcspeciesonly_6X75ind_poly
```



### Filtering for RAD2 dataset and subsequent analyses
### (fitlering across all populations; all individuals)

### bgzip then index file
```
bgzip MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20.recode.vcf
tabix MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20.recode.vcf.gz
```

### find sites with 75% having 6X cov USING ALL POPULATIONS
```
module load java/1.7
java -Xmx12G -jar GenomeAnalysisTK.jar \
-T CoveredByNSamplesSites \
-R pgl_assembly_v1.fasta \
--variant MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20.recode.vcf.gz \
-minCov 6 \
--percentageOfSamples .75 \
-out MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops \
-nt 12
```
### need to make tab deliminted
```
awk -F":" '$1=$1' OFS="\t" MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops > SNPs_hardfiltered_GQ20_6X75ind_allpops.txt
```

### filter for just these sites
```
vcftools --gzvcf MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20.recode.vcf.gz \
--positions SNPs_hardfiltered_GQ20_6X75ind_allpops.txt \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops
```
### calculate missingness by individual
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops.recode.vcf \
--missing-indv \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops
```
### calculate heterozygosity by individual
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops.recode.vcf \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops
```
### calculate mean depth by individual
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops.recode.vcf \
--depth \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops
```

### remove individuals with less than 14X MEAN coverage
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops.recode.vcf \
--keep List_14X_Meancoverage \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean
```


### Ouptut list of "RAD2" loci
```
output list of loci
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--site-depth \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean
```

### Assess sites out of HWE in the parental populations (P. glaucus)
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep Pg_species_list \
--hardy \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_Pg_HWE
```
### Assess sites out of HWE in the parental populations (P. canadensis)
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep Pc_species_list \
--hardy \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_Pc_HWE
```
### Combine output for both species populations
```
awk '{print > "PgandPc_HardGQ20_6X_31Xmean_HWE"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_Pg_HWE.hwe Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_Pc_HWE.hwe
```
### create a list of scaffolds to remove -- scaffolds with a p value < 0.01
```
awk '$6<0.01' PgandPc_HardGQ20_6X_31Xmean_HWE | sort | uniq > PgandPc_HardGQ20_6X_31Xmean_not_in_HWE
awk '{print $1" "$2}' PgandPc_HardGQ20_6X_31Xmean_not_in_HWE > PgandPc_HardGQ20_6X_31Xmean_sites_not_in_HWE
```
### Apply HWE filter
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--exclude-positions PgandPc_HardGQ20_6X_31Xmean_sites_not_in_HWE \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE
```


### Filter for sites that were identified as outliers using Bayescan and Arlequin - selected the 1st SNP for 1kb "windows"
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE.recode.vcf \
--positions Diversifying_outliers_1kb \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb
```
### print loci IDs for this filtered dataset
```
grep -v ^# Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf | awk '{print $1, $2}' > Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.loci
```


### Calculate allele freq for clinal loci for each population
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_1_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_1

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_2_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_2

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_3_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_3

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_4_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_4

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_5_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_5

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_6_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_6

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_7_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_7

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_HWE_outliers_1kb.recode.vcf \
--keep POP_8_NEW\
--freq \
--out SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_8
```

### print just alleles freq, remove header, split columns by semi colon, print just one of the allele frequencies (second column), transpose allele frequencies
```
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_1.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_1_SNPs_hard_6X_31Xmean_HWE_outliers_all
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_2.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_2_SNPs_hard_6X_31Xmean_HWE_outliers_all
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_3.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_3_SNPs_hard_6X_31Xmean_HWE_outliers_all
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_4.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_4_SNPs_hard_6X_31Xmean_HWE_outliers_all
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_5.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_5_SNPs_hard_6X_31Xmean_HWE_outliers_all
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_6.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_6_SNPs_hard_6X_31Xmean_HWE_outliers_all
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_7.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_7_SNPs_hard_6X_31Xmean_HWE_outliers_all
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_8.frq | cut -f 5,6 | tail -n +2 | awk '{ gsub(":", " ") } 1' | cut -f 2 | awk '{print $2}' | tr "\n" "\t" > AF_POP_8_SNPs_hard_6X_31Xmean_HWE_outliers_all
```

### Append them all together
```
awk '{print > "AF_All_hard_6X_31Xmean_HWE_outliers_all"}' AF_POP_1_SNPs_hard_6X_31Xmean_HWE_outliers_all AF_POP_2_SNPs_hard_6X_31Xmean_HWE_outliers_all AF_POP_3_SNPs_hard_6X_31Xmean_HWE_outliers_all AF_POP_4_SNPs_hard_6X_31Xmean_HWE_outliers_all AF_POP_5_SNPs_hard_6X_31Xmean_HWE_outliers_all AF_POP_6_SNPs_hard_6X_31Xmean_HWE_outliers_all AF_POP_7_SNPs_hard_6X_31Xmean_HWE_outliers_all AF_POP_8_SNPs_hard_6X_31Xmean_HWE_outliers_all
```

### Create vector of loci IDs (same loci used for all)
```
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_8.frq | cut -f 1,2 | tail -n +2 | tr "\t" "." | tr "\n" "\t" > AF_All_POPs_SNPs_hard_6X_31Xmean_HWE_outliers_all_Loci_IDs
```

### Determine number of individuals used to calc allele freq from above
```
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_1.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_1_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_2.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_2_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_3.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_3_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_4.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_4_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_5.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_5_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_6.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_6_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_7.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_7_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
cat SNPs_hard_6X_31Xmean_HWE_outliers_all_POP_8.frq | cut -f 4 | tail -n +2 | tr "\n" "\t" > AF_POP_8_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
```
### Now, need to append them all together
```
awk '{print > "All_POPs_SNPs_hard_6X_31Xmean_HWE_outliers_all_n"}' AF_POP_1_SNPs_hard_6X_31Xmean_HWE_outliers_all_n AF_POP_2_SNPs_hard_6X_31Xmean_HWE_outliers_all_n AF_POP_3_SNPs_hard_6X_31Xmean_HWE_outliers_all_n AF_POP_4_SNPs_hard_6X_31Xmean_HWE_outliers_all_n AF_POP_5_SNPs_hard_6X_31Xmean_HWE_outliers_all_n AF_POP_6_SNPs_hard_6X_31Xmean_HWE_outliers_all_n AF_POP_7_SNPs_hard_6X_31Xmean_HWE_outliers_all_n AF_POP_8_SNPs_hard_6X_31Xmean_HWE_outliers_all_n
```



### Scripts for processing Hzar output
### Remove spaces from file names  (need to be in BASH shell)
```
find . -name '* *' | while read file; do target=`echo "$file" | sed 's/ /_/g'`; echo "Renaming '$file' to '$target'"; mv "$file" "$target"; done;

mkdir MaxLL_var_params_for_selected_model
mkdir MaxLL_selected_model
mkdir MaxLL_params_for_selected_model
mkdir Selected_Models
#mkdir Rawdata_Plots
mkdir Tracemodels
mkdir CheckFit_related_files
mkdir Check_model_convergence
mkdir Fuzzyclines_for_selectedmodels
mkdir Comparing16models
mkdir AICtables
```
### Move files into appropriate folders
```
mv scaffold_* Results/
```

### move files with the same extension into a different directory
```
mv *AICc_table_for_all_models.txt Results/AICtables/
mv *tracemodel* Results/Tracemodels/
mv *rawdata.png Results/Rawdata_Plots/
mv *comparing16models.png Results/Comparing16models/
mv *fuzzycline_selectedmodel.png Results/Fuzzyclines_for_selectedmodels/
mv *__selected_model.txt Results/Selected_Models/
mv *_maxLL_selectedmodel.png Results/MaxLL_selected_model/
mv *_MaxLL_var_params_for_selected_model.txt Results/MaxLL_var_params_for_selected_model/
mv *_MaxLL_params_for_selected_model.txt Results/MaxLL_params_for_selected_model/
mv *_check* Results/CheckFit_related_files/
mv *_Check* Results/CheckFit_related_files/
```
### Merge all files (tables) in a single folder adding the file name to the first column in each row

### In AIC folder directory
```
grep "" *.txt > AICtable_for_6X_13Xmean_SNPs

#In Selected_Models directory
grep "" *.txt > Selected_6X_13Xmean_model

#In MaxLL_var_params_for_selected_model directory
grep "" *.txt > MaxLL_var_params_6X_13Xmean_selected_models

#In MaxLL_params_for_selected_model directory
grep "" *.txt > MaxLL_params_6X_13Xmean_selected_models
```








### Create dataset to use for determing genetic structure

### Filter for 1 kb SNPs (Rad2 dataset)
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--positions 6X_31Xmean_1kb \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_1kb
```
### print list of loci
```
grep -v ^# Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_1kb.recode.vcf | awk '{print $1, $2}' > Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_1kb.loci.txt
```
### convert vcf to structure format
```
module load java/1.7
java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar \
-inputfile Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_1kb.recode.vcf \
-inputformat VCF \
-outputfile Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_1kb.str \
-outputformat STRUCTURE \
-spid VCFtoSTR.spid
```

### Dataset for Introgress - loci with Fst > 0.9 (weir and cockerman) and 1 SNP per 1 kb

#### First, in excel make a list from the 31Xmean_1kb dataset and select only SNPs with Fst > 0.9
### Then, select only these sites
```
vcftools --vcf .recode.vcf \
--positions List_loci_for_introgress_6Xhard_1kb \
--recode \
--recode-INFO-all \
--out loci_for_introgress
```



### Calcualte r2 without using loop to be used to calucte LD decay

### Use sites > 100 bp apart
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--positions HardFiltered_SNPs_100bp \
--recode \
--recode-INFO-all \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp
```
### make plink file
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp.recode.vcf \
--plink \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink
```
### Add chrom numbers to map file in excel (a unique number to each scaffold)
#### label new file: "Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map"

### Convert plink file into 12coded plink file
```
plink-1.07-x86_64/plink --file Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom --recode12 --out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12
```
### make files with shorter names
```
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12.ped GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12.ped

cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_1_all_plinkfile12.map
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_2_all_plinkfile12.map
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_3_all_plinkfile12.map
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_4_all_plinkfile12.map
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_5_all_plinkfile12.map
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_6_all_plinkfile12.map
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_7_all_plinkfile12.map
cp Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom.map POP_8_all_plinkfile12.map
```

### Calculate r2 for each population
### clear all variables
```
rm POP_1_100bp_Iterative_r2_hard
rm POP_1_100bp_r2_mean_i
rm POP_1_100bp_r2.ld
rm POP_1_NEW_subset
rm Subsetted_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12_POP_1*
cp POP_1_all_plinkfile12.map Subsetted_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12_POP_1.map

#run loop to subsample from given population and calc r2 using plink, do this i times and record the mean each time, recording each mean of all scaffolds for each i

i="0"
while [ $i -lt 100 ]
do
        sort -R POP_1_NEW | head -n 7 > POP_1_NEW_subset

        grep -wFf POP_1_NEW_subset GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12.ped > Subsetted_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12_POP_1.ped
        pb-plink-src/plink --file Subsetted_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12_POP_1 --pb --r2 --ld-window-r2 0 --out POP_1_100bp_r2
        cat POP_1_100bp_r2.ld >> POP_1_100bp_Iterative_r2_hard
        i=$[$i+1]
done


rm POP_1_100bp_r2.ld
rm POP_1_NEW_subset
rm Subsetted_GQ20_6X75ind_allpops_31Xmean_100bp_plink_chrom_12_POP_1*

```

### process output for input into R
```
#create header
head -n 1 POP_1_100bp_Iterative_r2_hard > r2_dist_POPs_100bp_forR_header
```
### print columns to use and remove the multiple headers from each file
```
grep -v 'R2' POP_1_100bp_Iterative_r2_hard > r2_POP_1.1_100bp_forR
awk '{print > "r2_POP_1.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_1.1_100bp_forR

grep -v 'R2' POP_2_100bp_Iterative_r2_hard > r2_POP_2.1_100bp_forR
awk '{print > "r2_POP_2.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_2.1_100bp_forR

grep -v 'R2' POP_3_100bp_Iterative_r2_hard > r2_POP_3.1_100bp_forR
awk '{print > "r2_POP_3.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_3.1_100bp_forR

grep -v 'R2' POP_4_100bp_Iterative_r2_hard > r2_POP_4.1_100bp_forR
awk '{print > "r2_POP_4.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_4.1_100bp_forR

grep -v 'R2' POP_5_100bp_Iterative_r2_hard > r2_POP_5.1_100bp_forR
awk '{print > "r2_POP_5.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_5.1_100bp_forR

grep -v 'R2' POP_6_100bp_Iterative_r2_hard > r2_POP_6.1_100bp_forR
awk '{print > "r2_POP_6.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_6.1_100bp_forR

grep -v 'R2' POP_7_100bp_Iterative_r2_hard > r2_POP_7.1_100bp_forR
awk '{print > "r2_POP_7.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_7.1_100bp_forR

grep -v 'R2' POP_8_100bp_Iterative_r2_hard > r2_POP_8.1_100bp_forR
awk '{print > "r2_POP_8.1_100bp_forR.txt"}' r2_dist_POPs_100bp_forR_header r2_POP_8.1_100bp_forR
```










### Calculate Fis for each population
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_1_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP1

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_2_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP2

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_3_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP3

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_4_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP4

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_5_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP5

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_6_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP6

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_7_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP7

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_8_NEW \
--het \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP8
```
### comnbine them all together
```
cat Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP1.het Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP2.het Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP3.het Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP4.het Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP5.het Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP6.het Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP7.het Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP8.het > Het_6X_hard_31Xmean_all_pops
```





### Calculate pi for each population
```
vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_1_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP1

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_2_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP2

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_3_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP3

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_4_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP4

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_5_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP5

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_6_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP6

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_7_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP7

vcftools --vcf Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean.recode.vcf \
--keep POP_8_NEW \
--site-pi \
--out Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP8
```
###  now compile the output - add population info to each file and remove headers
```
awk '{print $3, "\t 1"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP1.sites.pi > pi_POP_1_hard
awk '{print $3, "\t 2"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP2.sites.pi > pi_POP_2_hard
awk '{print $3, "\t 3"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP3.sites.pi > pi_POP_3_hard
awk '{print $3, "\t 4"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP4.sites.pi > pi_POP_4_hard
awk '{print $3, "\t 5"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP5.sites.pi > pi_POP_5_hard
awk '{print $3, "\t 6"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP6.sites.pi > pi_POP_6_hard
awk '{print $3, "\t 7"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP7.sites.pi > pi_POP_7_hard
awk '{print $3, "\t 8"}' Variants_MiHiNew_NewOnly2_noindels_biallelic_less90miss_min05maf95_hardfiltered_GQ20_6X75ind_allpops_31Xmean_POP8.sites.pi > pi_POP_8_hard
```
### create header
```
head -n 1 pi_POP_1_hard | awk '{print $1, "\t Population"}' > pi_POPs_hard_header
```
### remove header from individual population files
```
awk '{if (NR!=1) {print}}' pi_POP_1_hard > pi_POP_1.1_hard
awk '{if (NR!=1) {print}}' pi_POP_2_hard > pi_POP_2.1_hard
awk '{if (NR!=1) {print}}' pi_POP_3_hard > pi_POP_3.1_hard
awk '{if (NR!=1) {print}}' pi_POP_4_hard > pi_POP_4.1_hard
awk '{if (NR!=1) {print}}' pi_POP_5_hard > pi_POP_5.1_hard
awk '{if (NR!=1) {print}}' pi_POP_6_hard > pi_POP_6.1_hard
awk '{if (NR!=1) {print}}' pi_POP_7_hard > pi_POP_7.1_hard
awk '{if (NR!=1) {print}}' pi_POP_8_hard > pi_POP_8.1_hard
```
### Combine all files together with a header
```
awk '{print > "pi_6X_hard_31Xmean_all_pops.txt"}' pi_POPs_hard_header pi_POP_1.1_hard pi_POP_2.1_hard pi_POP_3.1_hard pi_POP_4.1_hard pi_POP_5.1_hard pi_POP_6.1_hard pi_POP_7.1_hard pi_POP_8.1_hard
```