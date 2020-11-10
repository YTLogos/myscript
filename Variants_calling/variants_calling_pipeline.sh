#!/bin/bash

## This is the first part of the multi-sample variant calling pipeline, assuming that the data for each of your samples is illumina PE double-ended sequencing data fastq

set -e
set -u
set -o pipefail

## The path to the required software, the reference genome, and the file path.

bwa=/home/taoyan/biosoft/bwa/bwa
samtools=/home/taoyan/biosoft/samtools1.9/bin/samtools
reference=/database/reference/Bna_ZS11/ZS11_chr/ZS11_chr.fasta
gatk=/home/taoyan/biosoft/gatk-4.1.3.0/gatk
seq_data=/labdata/Bna_reseq_analysis/data/all_good_data
bam_data=/media/analysis/result/bam_data
vcf_data=/media/analysis/result/vcf_data

## In the first step, BWA establishes the reference genome index, enters the reference genome folder, and at the same time establishes the faidx index and dict file for subsequent GATK Calling SNP.

$bwa index $reference
$amtools faidx $reference
$gatk CreateSequenceDictionary \
-R $reference \
-O /database/reference/Bna_ZS11/ZS11_ZS11/ZS11_ZS11.dict && echo " **** dict done **** "


## The second step is the sequencing data quality control, all sample names are stored in a txt file, one line per sample name, sample_id.txt

## fastp is also OK 

#for i in `cat sample_id.txt`
#do
#$fastp -i $seq_data/${i}_1.fq.gz \
#-o $seq_data/${i}_good_1.fq.gz \
#-I $seq_data/${i}_2.fq.gz \
#-O $seq_data/${i}_good_2.fq.gz


for i in `cat sample_id.txt`
do
java -jar /home/taoyan/biosoft/Trimmomatic-0.30/trimmomatic-0.30.jar PE \ 
-threads 20 -phred33 $seq_data/${i}_1.fq.gz $seq_data/${i}_2.fq.gz \ 
$seq_data/${i}_good_1.fq.gz $seq_data/${i}_unpaired_1.fq.gz $seq_data/${i}_good_2.fq.gz $seq_data/${i}_unpaired_2.fq.gz \
ILLUMINACLIP:/home/taoyan/biosoft/Trimmomatic-0.30/adapters/TruSeq3-PE.fa:2:30:10 \ 
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
done

## The third step BWA MEM align, bwa mem is very effective for any reads longer than 40bp but less than 2000bp.
## To save time, here we use samtools to sort the generated files directly to generate the bam file
## Note that the number of threads and the amount of memory should not be set too large when sorting here, otherwise the error will be reported.

for i in `cat sample_id.txt`
do
$bwa mem -t 24 -R "@RG\tID:Brassica\tPL:illumina\tPU:reseq\tLB:Bna\tSM:${i}" \
$reference \
$seq_data/${i}_1.fq.gz $seq_data/${i}_2.fq.gz | $samtools view -Sb - > $bam_data/${i}.bam && \
echo "** Sample ${i} BWA MEM DONE **" && 
$samtools sort -@ 4 -m 8G -O bam -o $bam_data/${i}.sorted.bam $bam_data/${i}.bam && echo "** ${i}.bam sorted raw bamfile DONE **"

### Marked repeating sequences

$gatk MarkDuplicates \
   -I $bam_data/${i}.sorted.bam \
   -O $bam_data/${i}.sorted.markdup.bam \
   -M $bam_data/${i}.markdup_metrics.txt && echo "** ${i}.sorted.bam MarkDuplicates DONE **"

## Build Index for ${i}.sort.markdup.bam

$samtools index $bam_data/${i}.sorted.markdup.bam && echo "** ${i}.sorted.markdup.bam index DONE **"
rm -rf $bam_data/${i}.sorted.bam
rm -rf $bam_data/${i}.bam
done

## The fourth step is to call SNP.
### The gvcf of the output sample is implemented in two ways, with the same result but different speed.

##############################################################################################
### The first one outputs the full gvcf file of the sample directly and then Joint genotypeing

for i in `cat sample_id.txt`
do
$gatk --java-options "-Xmx16G -Djava.io.tmpdir=/database/tmp/" HaplotypeCaller \
--emit-ref-confidence GVCF \
-R $reference \
-I $bam_data/${i}.sorted.markdup.bam \
-O $vcf_data/${i}.HC.g.vcf.gz && echo " ****  ${i}.HC.g.vcf.gz done  **** "
done

## Joint genotypeing

sample_gvcfs=""
for i in `cat sample_id.txt`
do
sample_gvcfs=${sample_gvcfs}"-V $vcf_data/${i}.HC.g.vcf.gz \\" ## Note that this step needs to be changed to sample_gvcfs=${sample_gvcfs}"-V $vcf_data/${i}.HC.g.vcf.gz" if it's Centos OS, otherwise it won't be recognized.
done

$gatk --java-options "-Xmx16G -Djava.io.tmpdir=/database/tmp/" CombineGVCFs \
-R $reference \
${sample_gvcfs} \
-O $vcf_data/Bna.ZS11.HC.g.vcf.gz && echo " ***** Bna.ZS11.HC.g.vcf.gz done **** " && \
$gatk --java-options "-Xmx16G -Djava.io.tmpdir=/database/tmp/" GenotypeGVCFs \
-R $reference \
-V $vcf_data/Bna.ZS11.HC.g.vcf.gz \
-O $vcf_data/Bna.ZS11.HC.vcf.gz && echo " **** Bna.ZS11.HC.vcf.gz done **** "




#####################################################################################################
### The second outputs a vcf for each chromosome, and then merges the results from all the chromosomes to increase speed.

# First, separate the samples by chromosome and merge the gvcf of the corresponding chromosome for each sample.
# Then each will do Joint Calling on the chromosomes (GATK4 can only accept one gvcf input, so it needs to be merged first).

chrom=( A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 C01 C02 C03 C04 C05 C06 C07 C08 C09 )

for i in `cat sample_id.txt`
do
for chr in ${chrom[@]}
do
$gatk --java-options "-Xmx16G -Djava.io.tmpdir=/database/tmp/" HaplotypeCaller \
--emit-ref-confidence GVCF \
-R $reference \
-I $bam_data/${i}.sorted.markdup.bam \
-L $chr \
-O $vcf_data/${i}.HC.${chr}.g.vcf.gz && echo " ****  ${i}.HC.${chr}.g.vcf.gz done  **** " &
done 
done && wait

# Finally, Genotype results for each chromosome are merged, which requires that each sample must output gvcf by chromosome to increase speed

for chr in ${chrom[@]}
do
sample_gvcfs=""
for i in `cat sample_id.txt`
do
sample_gvcfs=${sample_gvcfs}"-V $vcf_data/${i}.HC.${chr}.g.vcf.gz \\"
done

$gatk --java-options "-Xmx16G -Djava.io.tmpdir=/database/tmp/" CombineGVCFs \
-R $refernce \
${sample_gvcfs} \
-O $vcf_data/Bna.ZS11.HC.${chr}.g.vcf.gz && echo " **** Bna.ZS11.HC.${chr}.g.vcf.gz done **** " && \
$gatk --java-options "-Xmx8G -Djava.io.tmpdir=/database/tmp/" GenotypeGVCFs \
-R $refernece \
-V $vcf_data/Bna.ZS11.HC.${chr}.g.vcf.gz \
-O $vcf_data/Bna.ZS11.Hc.${chr}.vcf.gz && echo " **** Bna.ZS11.HC.${chr}.vcf.gz done **** " &
done && wait

## MergeVcfs All chromosomes make up the final population VCF
merge_vcfs=""
for chr in ${chrom[@]}
do
merge_vcfs=${merge_vcfs}" -I $vcf_data/Bna.ZS11.Hc.${chr}.vcf.gz \\"
done && $gatk --java-options "-Xmx96G -Djava.io.tmpdir=/database/tmp/" MergeVcfs \
${merge_vcfs} \
-O $vcf_data/Bna.ZS11.HC.vcf.gz && echo " ****  MergeVcfs done  **** "

## Here's what separates SNPs and Indels

### SelectVariants, extract SNPs

$gatk --java-options "-Xmx120G -Djava.io.tmpdir=/database/tmp/" SelectVariants \
  -R $reference \
  -select-type SNP \
  -V $vcf_data/Bna.ZS11.HC.vcf.gz \
  -O $vcf_data/Bna.ZS11.HC.snp.raw.vcf.gz

### SNPs Filtering

$gatk --java-options "-Xms96g -Xmx120G -Djava.io.tmpdir=/database/tmp/" VariantFiltration \
  -R $reference \
  -V $vcf_data/Bna.ZS11.HC.snp.raw.vcf.gz \
  --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "Filter" \
  -O $vcf_data/Bna.ZS11.HC.snp.filter.vcf.gz
### Selecting SNPs with PASS tag
$gatk --java-options "-Xmx120G -Djava.io.tmpdir=/database/tmp/" SelectVariants \
-R $reference \
-V $vcf_data/Bna.ZS11.HC.snp.filter.vcf.gz \
--exclude-filtered \
-O $vcf_data/Bna.ZS11.snp.pass.vcf.gz


### SelectVariants, extract Indels
$gatk --java-options "-Xmx120G -Djava.io.tmpdir=/database/tmp/" SelectVariants \
  -R $reference \
  -select-type INDEL \
  -V $vcf_data/Bna.ZS11.HC.vcf.gz \
  -O $vcf_data/Bna.ZS11.HC.indel.raw.vcf.gz

### Indels Filtering
$gatk --java-options "-Xms96G -Xmx120G -Djava.io.tmpdir=/database/tmp/" VariantFiltration \
  -R $reference \
  -V $vcf_data/Bna.ZS11.HC.indel.raw.vcf.gz \
  --filter-expression "QUAL < 30.0 || QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
  --filter-name "Filter" \
  -O $vcf_data/Bna.ZS11.HC.indel.filter.vcf.gz

# Selecting INDELs with PASS tag
$gatk --java-options "-Xmx120G -Djava.io.tmpdir=/database/tmp/" SelectVariants \
-R $reference \
-V $vcf_data/Bna.ZS11.HC.indel.filter.vcf.gz \
--exclude-filtered \
-O $vcf_data/Bna.ZS11.indel.pass.vcf.gz

## Annotate SNPs

java -jar -Xms64G -Xmx96G "-Djava.io.tmpdir=/database/tmp/" /home/taoyan/biosoft/snpEff/snpEff.jar ZS11 $vcf_data/Bna.ZS11.HC.snp.pass.final.vcf.gz > $vcf_data/Bna.ZS11.HC.snp.pass.final.anno.vcf

/home/taoyan/biosoft/htslib1.9/bin/bgzip < $vcf_data/Bna.ZS11.HC.snp.pass.final.anno.vcf > $vcf_data/Bna.ZS11.HC.snp.pass.final.anno.vcf.gz

/home/taoyan/biosoft/htslib1.9/bin/tabix -p vcf $vcf_data/Bna.ZS11.HC.snp.pass.final.anno.vcf.gz

